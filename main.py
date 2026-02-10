"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D et Génétique
Base de données simulée de races ovines algériennes
Version avec critères de sélection mammaires et noms génériques
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
import json
import random
import math
import io
import base64
from PIL import Image, ImageDraw

# ============================================================================
# SECTION 2: CONFIGURATION STREAMLIT
# ============================================================================
st.set_page_config(
    page_title="Ovin Manager Pro - Races Algériennes",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# SECTION 3: CSS PERSONNALISÉ
# ============================================================================
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
    }
    .section-header {
        font-size: 2rem;
        color: #8B0000;
        margin-top: 2rem;
        margin-bottom: 1rem;
        padding-bottom: 10px;
        border-bottom: 3px solid #FF4500;
    }
    .race-card {
        background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(26,35,126,0.2);
    }
    .metric-card {
        background: linear-gradient(135deg, #FFF5F5 0%, #FFE4E1 100%);
        border-radius: 15px;
        padding: 15px;
        text-align: center;
        box-shadow: 0 5px 15px rgba(139,0,0,0.1);
        border-left: 5px solid #8B0000;
        margin: 5px;
    }
    .critere-card {
        background: linear-gradient(135deg, #2E7D32 0%, #4CAF50 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(46,125,50,0.2);
    }
    .mammelle-card {
        background: linear-gradient(135deg, #8B0000 0%, #FF4500 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(139,0,0,0.2);
    }
    .scanner-view {
        background: black;
        border-radius: 10px;
        padding: 10px;
        margin: 10px 0;
        text-align: center;
    }
    .gene-card {
        background: linear-gradient(135deg, #6a1b9a 0%, #8e24aa 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(138,27,154,0.2);
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# SECTION 4: STANDARDS DES RACES ALGÉRIENNES
# ============================================================================
STANDARDS_RACES = {
    'HAMRA': {
        'nom_complet': 'Hamra (Rousse)',
        'couleur': 'Rouge à marron',
        'origines': ['Sud Algérien', 'Sahara'],
        'caracteristiques': ['Robe rousse', 'Adaptée au désert', 'Bonne laitière'],
        'poids_adulte': {'femelle': (45, 65), 'male': (65, 90)},
        'mensurations': {
            'longueur_cm': (95, 125),
            'hauteur_cm': (65, 85),
            'tour_poitrine_cm': (95, 120),
            'largeur_bassin_cm': (35, 50)
        },
        'production_lait': (1.5, 3.5),
        'taux_mg': (6.0, 8.5),
        'prolificite': (1.2, 1.8),
        'criteres_mammaires': {
            'volume': 'Moyen à élevé',
            'trayons': '3-5 cm, bien orientés',
            'symetrie': 'Bonne',
            'aptitude_laitiere': 'Bonne'
        }
    },
    'OUDA': {
        'nom_complet': 'Ouled Djellal (Ouda)',
        'couleur': 'Blanche',
        'origines': ['Hauts Plateaux', 'Steppes'],
        'caracteristiques': ['Robe blanche', 'Queue grasse', 'Viande'],
        'poids_adulte': {'femelle': (50, 70), 'male': (70, 100)},
        'mensurations': {
            'longueur_cm': (100, 130),
            'hauteur_cm': (70, 90),
            'tour_poitrine_cm': (100, 130),
            'largeur_bassin_cm': (38, 55)
        },
        'production_lait': (1.0, 2.5),
        'taux_mg': (5.5, 7.5),
        'prolificite': (1.1, 1.5),
        'criteres_mammaires': {
            'volume': 'Grand',
            'trayons': '4-6 cm, légèrement divergents',
            'symetrie': 'Très bonne',
            'aptitude_laitiere': 'Excellente'
        }
    },
    'SIDAHOU': {
        'nom_complet': 'Sidahou',
        'couleur': 'Noire et blanche',
        'origines': ['Ouest Algérien'],
        'caracteristiques': ['Tête noire', 'Résistante', 'Mixte'],
        'poids_adulte': {'femelle': (40, 60), 'male': (60, 85)},
        'mensurations': {
            'longueur_cm': (90, 120),
            'hauteur_cm': (60, 80),
            'tour_poitrine_cm': (90, 115),
            'largeur_bassin_cm': (34, 48)
        },
        'production_lait': (1.2, 2.8),
        'taux_mg': (6.2, 8.0),
        'prolificite': (1.3, 1.7),
        'criteres_mammaires': {
            'volume': 'Moyen',
            'trayons': '3-5 cm, bien insérés',
            'symetrie': 'Bonne',
            'aptitude_laitiere': 'Moyenne'
        }
    },
    'BERBERE': {
        'nom_complet': 'Brebis Berbère',
        'couleur': 'Variée',
        'origines': ['Kabylie', 'Aurès'],
        'caracteristiques': ['Rustique', 'Petite taille', 'Adaptée montagne'],
        'poids_adulte': {'femelle': (35, 50), 'male': (50, 70)},
        'mensurations': {
            'longueur_cm': (80, 110),
            'hauteur_cm': (55, 75),
            'tour_poitrine_cm': (85, 110),
            'largeur_bassin_cm': (30, 45)
        },
        'production_lait': (0.8, 2.0),
        'taux_mg': (6.5, 9.0),
        'prolificite': (1.0, 1.4),
        'criteres_mammaires': {
            'volume': 'Petit à moyen',
            'trayons': '2-4 cm, bien insérés',
            'symetrie': 'Moyenne',
            'aptitude_laitiere': 'Adaptée'
        }
    },
    'CROISE': {
        'nom_complet': 'Croisement',
        'couleur': 'Variable',
        'origines': ['Multiple'],
        'caracteristiques': ['Vigueur hybride', 'Adaptabilité'],
        'poids_adulte': {'femelle': (40, 70), 'male': (60, 95)},
        'mensurations': {
            'longueur_cm': (85, 125),
            'hauteur_cm': (60, 85),
            'tour_poitrine_cm': (90, 125),
            'largeur_bassin_cm': (32, 52)
        },
        'production_lait': (1.0, 3.0),
        'taux_mg': (5.5, 8.5),
        'prolificite': (1.2, 1.8),
        'criteres_mammaires': {
            'volume': 'Variable',
            'trayons': '3-6 cm, orientation variable',
            'symetrie': 'Variable',
            'aptitude_laitiere': 'Variable'
        }
    },
    'INCONNU': {
        'nom_complet': 'Race non identifiée',
        'couleur': 'Indéterminée',
        'origines': ['Inconnue'],
        'caracteristiques': ['À caractériser'],
        'poids_adulte': {'femelle': (30, 60), 'male': (50, 80)},
        'mensurations': {
            'longueur_cm': (80, 120),
            'hauteur_cm': (55, 80),
            'tour_poitrine_cm': (85, 120),
            'largeur_bassin_cm': (30, 50)
        },
        'production_lait': (0.5, 2.5),
        'taux_mg': (5.0, 8.0),
        'prolificite': (1.0, 1.6),
        'criteres_mammaires': {
            'volume': 'À évaluer',
            'trayons': 'À mesurer',
            'symetrie': 'À évaluer',
            'aptitude_laitiere': 'À déterminer'
        }
    }
}

# ============================================================================
# SECTION 5: FONCTIONS STATISTIQUES (sans scipy)
# ============================================================================
def skewness(data):
    """Calcule le coefficient d'asymétrie de Pearson"""
    if len(data) < 3:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 3)

def kurtosis(data):
    """Calcule le coefficient d'aplatissement"""
    if len(data) < 4:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 4) - 3

# ============================================================================
# SECTION 6: BASE DE DONNÉES
# ============================================================================
def creer_base_races():
    """Crée une base de données avec races algériennes"""
    conn = sqlite3.connect('ovin_algerien.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table brebis avec critères mammaires
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE NOT NULL,
            nom TEXT,
            race TEXT CHECK(race IN ('HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU')),
            sous_race TEXT,
            sexe TEXT CHECK(sexe IN ('F', 'M')),
            date_naissance DATE,
            age_mois INTEGER,
            poids FLOAT,
            score_condition INTEGER CHECK(score_condition BETWEEN 1 AND 5),
            
            -- Caractères morphologiques
            couleur_robe TEXT,
            intensite_couleur INTEGER CHECK(intensite_couleur BETWEEN 1 AND 10),
            cornes BOOLEAN,
            taille_cornes_cm FLOAT,
            forme_cornes TEXT,
            type_laine TEXT,
            qualite_laine INTEGER CHECK(qualite_laine BETWEEN 1 AND 10),
            
            -- Mensurations
            longueur_corps_cm FLOAT,
            hauteur_garrot_cm FLOAT,
            largeur_bassin_cm FLOAT,
            tour_poitrine_cm FLOAT,
            circonference_tete_cm FLOAT,
            longueur_oreille_cm FLOAT,
            
            -- Critères mammaires (nouveau)
            volume_mammaire INTEGER CHECK(volume_mammaire BETWEEN 1 AND 5),
            symetrie_mammaire INTEGER CHECK(symetrie_mammaire BETWEEN 1 AND 5),
            insertion_trayons INTEGER CHECK(insertion_trayons BETWEEN 1 AND 5),
            longueur_trayons_cm FLOAT,
            orientation_trayons TEXT,
            
            -- Caractères qualitatifs
            temperement TEXT CHECK(temperement IN ('calme', 'nervieux', 'intermediaire')),
            aptitude TEXT CHECK(aptitude IN ('lait', 'viande', 'mixte', 'laine')),
            score_conformation FLOAT,
            aptitudes TEXT,
            notes TEXT,
            
            -- Génétique
            mere_id TEXT,
            pere_id TEXT,
            coefficient_consanguinite FLOAT DEFAULT 0.0,
            statut TEXT DEFAULT 'active',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Table production laitière
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS production_lait (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_mesure DATE,
            quantite_litre FLOAT,
            taux_matiere_grasse FLOAT,
            taux_proteine FLOAT,
            cellules_somatiques INTEGER,
            lactose FLOAT,
            ph FLOAT,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Table scanner 3D
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS scans_3d (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_scan DATE,
            mode_scan TEXT,
            points_3d_json TEXT,
            mesures_json TEXT,
            volume_estime FLOAT,
            surface_estimee FLOAT,
            qualite_scan INTEGER,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Table génotypage
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS genotypage (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            marqueur TEXT,
            chromosome TEXT,
            position INTEGER,
            allele1 TEXT,
            allele2 TEXT,
            genotype TEXT,
            frequence_allelique FLOAT,
            effet_additif FLOAT,
            effet_dominant FLOAT,
            r2 FLOAT,
            p_value FLOAT,
            gene_associe TEXT,
            trait_associe TEXT,
            date_analyse DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Vérifier si la base est vide et la peupler
    cursor.execute("SELECT COUNT(*) FROM brebis")
    count = cursor.fetchone()[0]
    
    if count == 0:
        peupler_base_races(cursor, conn)
    
    conn.commit()
    return conn

def peupler_base_races(cursor, conn):
    """Peuple la base avec des races algériennes - AVEC NOMS GÉNÉRIQUES"""
    races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU']
    brebis_data = []
    
    for i in range(1, 51):
        race = random.choice(races)
        sexe = random.choices(['F', 'M'], weights=[0.7, 0.3])[0]
        
        # Générer identifiant et nom générique (SANS NOMS ARABES/FRANÇAIS)
        identifiant = f"{race[:3]}-{sexe}-2023-{i:03d}"
        
        # Noms génériques : F pour femelle, M pour mâle + race + numéro
        if sexe == 'F':
            nom = f"F{race[:3]}{i:03d}"  # Exemple: FHAM001, FOUDA002
        else:
            nom = f"M{race[:3]}{i:03d}"  # Exemple: MHAM001, MOUDA002
        
        age_mois = random.randint(12, 84)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        
        # Poids selon race et sexe
        poids_min, poids_max = STANDARDS_RACES[race]['poids_adulte'][sexe.lower()]
        poids = random.uniform(poids_min, poids_max)
        
        # Score de condition
        score_condition = random.randint(2, 4)
        
        # Couleur selon race
        couleurs = {
            'HAMRA': ['Rousse', 'Rousse foncée', 'Marron'],
            'OUDA': ['Blanche', 'Crème', 'Blanc cassé'],
            'SIDAHOU': ['Noire et blanche', 'Pie noire', 'Tête noire'],
            'BERBERE': ['Noire', 'Brune', 'Grise', 'Pie'],
            'CROISE': ['Variable', 'Panachée', 'Mélangée'],
            'INCONNU': ['Indéterminée']
        }
        couleur_robe = random.choice(couleurs[race])
        
        # Critères mammaires (seulement pour femelles)
        if sexe == 'F':
            volume_mammaire = random.randint(2, 5)
            symetrie_mammaire = random.randint(2, 5)
            insertion_trayons = random.randint(2, 5)
            longueur_trayons = random.uniform(3.0, 6.0)
            orientation_trayons = random.choice(['parallele', 'leger_divergent', 'divergent'])
        else:
            volume_mammaire = None
            symetrie_mammaire = None
            insertion_trayons = None
            longueur_trayons = None
            orientation_trayons = None
        
        # Score conformation calculé
        score_conformation = random.uniform(5.0, 9.0)
        
        # Mensurations
        longueur_corps = random.uniform(*STANDARDS_RACES[race]['mensurations']['longueur_cm'])
        hauteur_garrot = random.uniform(*STANDARDS_RACES[race]['mensurations']['hauteur_cm'])
        largeur_bassin = random.uniform(*STANDARDS_RACES[race]['mensurations']['largeur_bassin_cm'])
        tour_poitrine = random.uniform(*STANDARDS_RACES[race]['mensurations']['tour_poitrine_cm'])
        
        brebis_data.append((
            identifiant, nom, race, '', sexe, date_naissance.isoformat(), 
            age_mois, poids, score_condition, couleur_robe, 
            random.randint(5, 10), random.choice([True, False]), 
            random.uniform(0, 60), '', random.choice(['fine', 'semi-fine', 'grossière']), 
            random.randint(3, 9),
            longueur_corps, hauteur_garrot, largeur_bassin, tour_poitrine,
            random.uniform(45, 65), random.uniform(12, 18),
            volume_mammaire, symetrie_mammaire, insertion_trayons,
            longueur_trayons, orientation_trayons,
            random.choice(['calme', 'nervieux', 'intermediaire']),
            random.choice(['lait', 'viande', 'mixte', 'laine']),
            score_conformation, '',
            f"Brebis {race} - Élevage algérien",
            None, None, random.uniform(0.0, 0.15), 'active'
        ))
    
    # Insérer les brebis
    cursor.executemany('''
        INSERT INTO brebis (
            identifiant, nom, race, sous_race, sexe, date_naissance, age_mois, 
            poids, score_condition, couleur_robe, intensite_couleur, cornes, 
            taille_cornes_cm, forme_cornes, type_laine, qualite_laine,
            longueur_corps_cm, hauteur_garrot_cm, largeur_bassin_cm, 
            tour_poitrine_cm, circonference_tete_cm, longueur_oreille_cm,
            volume_mammaire, symetrie_mammaire, insertion_trayons,
            longueur_trayons_cm, orientation_trayons,
            temperement, aptitude, score_conformation, aptitudes, notes, 
            mere_id, pere_id, coefficient_consanguinite, statut
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', brebis_data)
    
    conn.commit()

# Initialiser la base
conn = creer_base_races()

# ============================================================================
# SECTION 7: MODULE SCANNER 3D
# ============================================================================
class Scanner3D:
    """Simulateur de scanner 3D pour ovins"""
    
    @staticmethod
    def generer_photo_simulee(brebis_info):
        """Génère une photo simulée selon la race"""
        width, height = 400, 300
        image = Image.new('RGB', (width, height), color='white')
        draw = ImageDraw.Draw(image)
        
        # Couleur selon la race
        couleurs = {
            'HAMRA': (139, 0, 0),      # Rouge foncé
            'OUDA': (255, 255, 255),   # Blanc
            'SIDAHOU': (50, 50, 50),   # Noir
            'BERBERE': (165, 42, 42),  # Brun
            'CROISE': (160, 120, 80),  # Marron
            'INCONNU': (200, 200, 200) # Gris
        }
        
        race = brebis_info.get('race', 'INCONNU')
        corps_color = couleurs.get(race, (200, 200, 200))
        
        # Dessiner le corps
        draw.ellipse([100, 80, 300, 200], fill=corps_color, outline='black', width=2)
        
        # Tête
        draw.ellipse([280, 100, 350, 160], fill=corps_color, outline='black', width=2)
        
        # Pattes
        for x in [130, 170, 230, 270]:
            draw.rectangle([x, 200, x+20, 280], fill='black')
        
        # Cornes si présentes
        if brebis_info.get('cornes', False):
            draw.arc([300, 70, 350, 120], start=0, end=180, fill='gray', width=5)
            draw.arc([320, 70, 370, 120], start=0, end=180, fill='gray', width=5)
        
        # Informations
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {race}", fill='black')
        draw.text((10, 50), f"Poids: {brebis_info.get('poids', 0):.1f} kg", fill='black')
        
        return image
    
    @staticmethod
    def simuler_scan_3d(brebis_info):
        """Simule un scan 3D réaliste"""
        np.random.seed(hash(brebis_info.get('identifiant', '')) % 10000)
        
        n_points = 500
        points = []
        
        poids = brebis_info.get('poids', 50)
        race = brebis_info.get('race', 'INCONNU')
        
        # Rayons approximatifs
        rx = 0.6 * poids**0.33  # Largeur
        ry = 1.2 * poids**0.33  # Longueur
        rz = 0.8 * poids**0.33  # Hauteur
        
        for _ in range(n_points):
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
            
            x = rx * np.sin(phi) * np.cos(theta) + np.random.normal(0, rx*0.05)
            y = ry * np.sin(phi) * np.sin(theta) + np.random.normal(0, ry*0.05)
            z = rz * np.cos(phi) + np.random.normal(0, rz*0.05)
            
            # Intensité selon la position
            if z > rz * 0.5:
                intensity = np.random.uniform(100, 150)
            elif abs(x) > rx * 0.7:
                intensity = np.random.uniform(150, 200)
            else:
                intensity = np.random.uniform(200, 255)
            
            points.append({
                'x': float(x),
                'y': float(y),
                'z': float(z),
                'intensity': int(intensity)
            })
        
        return points

# ============================================================================
# SECTION 8: MODULE GÉNÉTIQUE
# ============================================================================
class ModuleGenetique:
    """Module d'analyse génétique"""
    
    @staticmethod
    def generer_genotype(brebis_id, race):
        """Génère un génotype simulé"""
        genotypes = []
        
        for i in range(10):
            marqueur = f"SNP{i+1:03d}"
            chromosome = str(random.randint(1, 26))
            position = random.randint(1000000, 90000000)
            allele1 = random.choice(['A', 'C', 'G', 'T'])
            allele2 = random.choice(['A', 'C', 'G', 'T'])
            genotype = allele1 + allele2
            
            genotypes.append((
                brebis_id, marqueur, chromosome, position, allele1, allele2,
                genotype, random.uniform(0.1, 0.9), random.uniform(-0.5, 0.5),
                random.uniform(-0.3, 0.3), random.uniform(0.1, 0.3),
                random.uniform(0.001, 0.05), f"GENE_{marqueur}",
                random.choice(['poids', 'production_lait', 'couleur', 'resistance']),
                date.today().isoformat()
            ))
        
        return genotypes
    
    @staticmethod
    def calculer_diversite_genetique(genotypes):
        """Calcule la diversité génétique - VERSION SÉCURISÉE"""
        if not genotypes:
            return {}
        
        # Convertir en liste de dictionnaires de manière sécurisée
        data = []
        for geno in genotypes:
            if len(geno) >= 8:  # Minimum pour les calculs de base
                data.append({
                    'marqueur': geno[1] if len(geno) > 1 else '',
                    'allele1': geno[4] if len(geno) > 4 else '',
                    'allele2': geno[5] if len(geno) > 5 else '',
                    'freq_allelique': float(geno[7]) if len(geno) > 7 else 0.5
                })
        
        if not data:
            return {}
        
        df = pd.DataFrame(data)
        
        # Calculs de diversité
        if 'allele1' in df.columns and 'allele2' in df.columns:
            heterozygotes = df[df['allele1'] != df['allele2']]
            ho = len(heterozygotes) / len(df) if len(df) > 0 else 0
        else:
            ho = 0
        
        if 'freq_allelique' in df.columns:
            he = 1 - (df['freq_allelique']**2).mean()
        else:
            he = 0
        
        fis = 1 - (ho / he) if he > 0 else 0
        
        return {
            'heterozygosite_observee': round(ho, 4),
            'heterozygosite_attendue': round(he, 4),
            'fis': round(fis, 4),
            'nombre_snps': len(df['marqueur'].unique()) if 'marqueur' in df.columns else 0
        }

# ============================================================================
# SECTION 9: PAGE ACCUEIL
# ============================================================================
def page_accueil():
    """Page d'accueil avec vue d'ensemble"""
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - RACES ALGÉRIENNES</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique des races ovines algériennes**")
    
    # Métriques principales
    cursor = conn.cursor()
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-card'>
            <h3>🐑 TOTAL BREBIS</h3>
            <h2>{total}</h2>
            <p>Races algériennes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        cursor.execute("SELECT COUNT(DISTINCT race) FROM brebis")
        races = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-card'>
            <h3>🏷️ RACES</h3>
            <h2>{races}</h2>
            <p>Différentes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        cursor.execute("SELECT AVG(poids) FROM brebis WHERE sexe = 'F'")
        poids_f = cursor.fetchone()[0] or 0
        st.markdown(f"""
        <div class='metric-card'>
            <h3>♀️ POIDS MOYEN</h3>
            <h2>{poids_f:.1f} kg</h2>
            <p>Femelles</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        cursor.execute("SELECT COUNT(*) FROM scans_3d")
        scans = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-card'>
            <h3>📐 SCANS 3D</h3>
            <h2>{scans}</h2>
            <p>Réalisés</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Distribution des races
    st.markdown("### 📊 DISTRIBUTION DES RACES")
    
    cursor.execute("""
        SELECT race, COUNT(*) as count,
               AVG(poids) as poids_moyen,
               AVG(age_mois) as age_moyen
        FROM brebis
        GROUP BY race
        ORDER BY count DESC
    """)
    
    races_data = cursor.fetchall()
    
    if races_data:
        df_races = pd.DataFrame(races_data, columns=['Race', 'Nombre', 'Poids moyen', 'Âge moyen'])
        
        col_race1, col_race2 = st.columns([2, 1])
        
        with col_race1:
            fig = px.pie(df_races, values='Nombre', names='Race',
                        title="Répartition des races dans le troupeau",
                        hole=0.4,
                        color_discrete_sequence=px.colors.sequential.RdBu)
            st.plotly_chart(fig, use_container_width=True)
        
        with col_race2:
            st.markdown("### 📈 CARACTÉRISTIQUES PAR RACE")
            for _, row in df_races.iterrows():
                race_info = STANDARDS_RACES.get(row['Race'], {})
                st.markdown(f"""
                <div class='race-card'>
                    <h4>{race_info.get('nom_complet', row['Race'])}</h4>
                    <p><strong>{row['Nombre']}</strong> animaux</p>
                    <p>Poids moyen: <strong>{row['Poids moyen']:.1f} kg</strong></p>
                    <p>Âge moyen: <strong>{row['Âge moyen']:.0f} mois</strong></p>
                </div>
                """, unsafe_allow_html=True)

# ============================================================================
# SECTION 10: PAGE SCANNER 3D
# ============================================================================
def page_scanner_3d():
    """Page du scanner 3D avec saisie manuelle"""
    st.markdown('<h2 class="section-header">📐 SCANNER 3D & SAISIE MANUELLE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🎯 SCANNER 3D", "📝 SAISIE MANUELLE"])
    
    with tab1:
        # Sélection de la brebis
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY nom")
        brebis_list = cursor.fetchall()
        
        if not brebis_list:
            st.warning("Aucune brebis dans la base de données")
            return
        
        col_sel1, col_sel2 = st.columns([2, 1])
        
        with col_sel1:
            brebis_option = st.selectbox(
                "SÉLECTIONNEZ UNE BREBIS À SCANNER:",
                [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list],
                key="scanner_selection"
            )
        
        with col_sel2:
            mode_scan = st.selectbox(
                "MODE DE SCAN:",
                ["Laser haute précision", "Photogrammétrie", "Scanner portable", "Simulation"],
                key="mode_scan"
            )
        
        if brebis_option:
            try:
                brebis_id = int(brebis_option.split('(')[1].split(')')[0].split('-')[-1])
            except:
                st.error("Erreur dans le format de l'identifiant")
                return
            
            cursor.execute("SELECT * FROM brebis WHERE id = ?", (brebis_id,))
            brebis_info = cursor.fetchone()
            columns = [desc[0] for desc in cursor.description]
            brebis_dict = dict(zip(columns, brebis_info))
            
            # Onglets du scanner
            scan_tabs = st.tabs(["📸 PHOTO", "🎯 SCAN 3D", "📏 MESURES"])
            
            with scan_tabs[0]:
                photo = Scanner3D.generer_photo_simulee(brebis_dict)
                
                col_photo1, col_photo2 = st.columns([2, 1])
                
                with col_photo1:
                    st.image(photo, caption=f"Photo simulée - {brebis_dict['nom']}", use_column_width=True)
                    
                    col_btn1, col_btn2 = st.columns(2)
                    with col_btn1:
                        if st.button("📸 Prendre photo", type="primary", key="photo_btn"):
                            st.success("Photo prise et sauvegardée!")
                    
                    with col_btn2:
                        if st.button("🔄 Régénérer", key="regenerate_btn"):
                            st.rerun()
                
                with col_photo2:
                    st.markdown(f"""
                    <div class='race-card'>
                        <h4>📷 INFORMATIONS</h4>
                        <p><strong>Race:</strong> {brebis_dict['race']}</p>
                        <p><strong>ID:</strong> {brebis_dict['identifiant']}</p>
                        <p><strong>Sexe:</strong> {brebis_dict['sexe']}</p>
                        <p><strong>Âge:</strong> {brebis_dict['age_mois']} mois</p>
                        <p><strong>Poids:</strong> {brebis_dict['poids']:.1f} kg</p>
                    </div>
                    """, unsafe_allow_html=True)
            
            with scan_tabs[1]:
                st.markdown("### 🎯 SCAN 3D EN TEMPS RÉEL")
                
                scan_progress = st.slider("Progression du scan:", 0, 100, 50, key="scan_progress")
                
                if scan_progress > 30:
                    points = Scanner3D.simuler_scan_3d(brebis_dict)
                    
                    st.markdown("**Points 3D capturés:**")
                    df_points = pd.DataFrame(points[:10])
                    st.dataframe(df_points[['x', 'y', 'z', 'intensity']])
                    
                    if st.button("🚀 Démarrer scan complet", type="primary", key="scan_btn"):
                        with st.spinner("Scan en cours... Veuillez patienter"):
                            for i in range(1, 101):
                                st.progress(i/100)
                            
                            cursor.execute('''
                                INSERT INTO scans_3d (brebis_id, date_scan, mode_scan, points_3d_json, qualite_scan, notes)
                                VALUES (?, ?, ?, ?, ?, ?)
                            ''', (
                                brebis_id,
                                date.today().isoformat(),
                                mode_scan,
                                json.dumps(points[:200]),
                                85,
                                f"Scan {mode_scan} - {brebis_dict['nom']}"
                            ))
                            conn.commit()
                            st.success("✅ Scan 3D terminé et sauvegardé!")
    
    with tab2:
        st.markdown("### 📝 SAISIE MANUELLE DES MESURES")
        
        with st.form("saisie_manuelle"):
            col_saisie1, col_saisie2 = st.columns(2)
            
            with col_saisie1:
                identifiant = st.text_input("Identifiant de l'animal")
                date_mesure = st.date_input("Date de mesure", value=date.today())
                operateur = st.text_input("Opérateur")
                
                st.markdown("#### 📏 MENSURATIONS (cm)")
                longueur = st.number_input("Longueur corps", 0.0, 200.0, 100.0, 0.1)
                hauteur = st.number_input("Hauteur garrot", 0.0, 150.0, 70.0, 0.1)
                largeur = st.number_input("Largeur bassin", 0.0, 100.0, 40.0, 0.1)
                tour_poitrine = st.number_input("Tour de poitrine", 0.0, 200.0, 100.0, 0.1)
            
            with col_saisie2:
                poids = st.number_input("Poids (kg)", 0.0, 200.0, 50.0, 0.1)
                score_condition = st.slider("Score de condition", 1, 5, 3)
                
                st.markdown("#### 🎨 CARACTÈRES QUALITATIFS")
                couleur_robe = st.text_input("Couleur de la robe")
                etat_corporel = st.select_slider("État corporel", 
                                                ['Maigre', 'Normal', 'Gras'])
                temperement = st.selectbox("Tempérament", 
                                         ['Calme', 'Nervieux', 'Intermediaire', 'Agité'])
                notes = st.text_area("Notes complémentaires")
            
            if st.form_submit_button("💾 Enregistrer les mesures", type="primary"):
                # Calculer des indices
                indice_corporel = (longueur * hauteur * largeur) ** (1/3)
                ratio_conformation = longueur / hauteur
                
                st.success(f"✅ Mesures enregistrées pour {identifiant}")
                st.info(f"""
                **Résumé:**
                - Indice corporel: {indice_corporel:.1f}
                - Ratio conformation: {ratio_conformation:.2f}
                - Poids: {poids} kg
                - Score: {score_condition}/5
                """)

# ============================================================================
# SECTION 11: PAGE GESTION
# ============================================================================
def page_gestion():
    """Page de gestion du troupeau"""
    st.markdown('<h2 class="section-header">📊 GESTION DU TROUPEAU</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🐑 LISTE", "📈 STATISTIQUES", "🔍 RECHERCHE", "📤 EXPORT"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT identifiant, nom, race, sexe, age_mois, poids, 
                   score_condition, couleur_robe, statut
            FROM brebis
            ORDER BY race, identifiant
        """)
        
        brebis_data = cursor.fetchall()
        columns = ['ID', 'Nom', 'Race', 'Sexe', 'Âge', 'Poids', 'Score', 'Couleur', 'Statut']
        
        df = pd.DataFrame(brebis_data, columns=columns)
        
        # Filtres
        col_filtre1, col_filtre2, col_filtre3 = st.columns(3)
        
        with col_filtre1:
            races_select = st.multiselect("Filtrer par race", df['Race'].unique(), df['Race'].unique())
        
        with col_filtre2:
            sexe_select = st.multiselect("Filtrer par sexe", df['Sexe'].unique(), df['Sexe'].unique())
        
        with col_filtre3:
            recherche = st.text_input("Recherche texte")
        
        # Appliquer filtres
        if races_select:
            df = df[df['Race'].isin(races_select)]
        if sexe_select:
            df = df[df['Sexe'].isin(sexe_select)]
        if recherche:
            df = df[df.apply(lambda row: row.astype(str).str.contains(recherche, case=False).any(), axis=1)]
        
        st.dataframe(df, use_container_width=True, height=500)
        
        st.metric("Brebis affichées", len(df), f"Sur {len(brebis_data)} au total")
    
    with tab2:
        st.markdown("### 📊 STATISTIQUES PAR RACE")
        
        cursor.execute("""
            SELECT race,
                   COUNT(*) as nombre,
                   ROUND(AVG(poids), 1) as poids_moyen,
                   ROUND(AVG(age_mois), 0) as age_moyen,
                   ROUND(AVG(score_condition), 1) as score_moyen,
                   ROUND(AVG(longueur_corps_cm), 1) as longueur_moyenne,
                   ROUND(AVG(hauteur_garrot_cm), 1) as hauteur_moyenne
            FROM brebis
            GROUP BY race
            ORDER BY nombre DESC
        """)
        
        stats_data = cursor.fetchall()
        
        if stats_data:
            df_stats = pd.DataFrame(stats_data, 
                                  columns=['Race', 'Nombre', 'Poids moyen', 'Âge moyen', 
                                           'Score moyen', 'Longueur moyenne', 'Hauteur moyenne'])
            
            # Graphiques
            col_stat1, col_stat2 = st.columns(2)
            
            with col_stat1:
                fig = px.bar(df_stats, x='Race', y='Poids moyen',
                            title="Poids moyen par race",
                            color='Nombre',
                            color_continuous_scale='Reds')
                st.plotly_chart(fig, use_container_width=True)
            
            with col_stat2:
                fig = px.scatter(df_stats, x='Longueur moyenne', y='Hauteur moyenne',
                                size='Nombre', color='Race',
                                title="Relation longueur/hauteur par race",
                                hover_data=['Race', 'Nombre', 'Poids moyen'])
                st.plotly_chart(fig, use_container_width=True)
            
            # Tableau détaillé
            st.markdown("### 📋 TABLEAU DÉTAILLÉ")
            st.dataframe(df_stats.style.background_gradient(subset=['Poids moyen'], cmap='Reds'))
    
    with tab3:
        st.markdown("### 🔍 RECHERCHE AVANCÉE")
        
        with st.form("recherche_avancee"):
            col_rech1, col_rech2 = st.columns(2)
            
            with col_rech1:
                min_poids = st.number_input("Poids minimum (kg)", 0, 200, 30)
                max_poids = st.number_input("Poids maximum (kg)", 0, 200, 100)
                min_score = st.slider("Score condition minimum", 1, 5, 2)
            
            with col_rech2:
                races = list(STANDARDS_RACES.keys())
                races_select = st.multiselect("Races", races, races)
                avec_cornes = st.selectbox("Cornes", ["Tous", "Avec", "Sans"])
                temperement = st.selectbox("Tempérament", ["Tous", "calme", "nervieux", "intermediaire"])
            
            if st.form_submit_button("🔍 Rechercher"):
                query = "SELECT * FROM brebis WHERE 1=1"
                params = []
                
                query += " AND poids BETWEEN ? AND ?"
                params.extend([min_poids, max_poids])
                
                query += " AND score_condition >= ?"
                params.append(min_score)
                
                if races_select:
                    placeholders = ','.join(['?'] * len(races_select))
                    query += f" AND race IN ({placeholders})"
                    params.extend(races_select)
                
                if avec_cornes == "Avec":
                    query += " AND cornes = 1"
                elif avec_cornes == "Sans":
                    query += " AND cornes = 0"
                
                if temperement != "Tous":
                    query += " AND temperement = ?"
                    params.append(temperement)
                
                cursor.execute(query, params)
                resultats = cursor.fetchall()
                
                if resultats:
                    df_result = pd.DataFrame(resultats, columns=[desc[0] for desc in cursor.description])
                    st.success(f"🔍 {len(resultats)} résultats trouvés")
                    st.dataframe(df_result[['identifiant', 'nom', 'race', 'sexe', 'poids', 'score_condition', 'couleur_robe']])
                else:
                    st.warning("Aucun résultat trouvé")
    
    with tab4:
        st.markdown("### 📤 EXPORTATION DES DONNÉES")
        
        export_type = st.selectbox("Type de données à exporter", 
                                  ["Troupeau complet", "Données morphologiques", "Production laitière", "Scans 3D"])
        
        export_format = st.selectbox("Format", ["CSV", "JSON"])
        
        if st.button("📥 Générer l'export", type="primary"):
            cursor.execute("SELECT * FROM brebis")
            data = cursor.fetchall()
            columns = [desc[0] for desc in cursor.description]
            df_export = pd.DataFrame(data, columns=columns)
            
            if export_format == "CSV":
                csv = df_export.to_csv(index=False, encoding='utf-8-sig')
                st.download_button(
                    label="📥 Télécharger CSV",
                    data=csv,
                    file_name=f"troupeau_{date.today()}.csv",
                    mime="text/csv"
                )
            elif export_format == "JSON":
                json_data = df_export.to_json(orient='records', indent=2, force_ascii=False)
                st.download_button(
                    label="📥 Télécharger JSON",
                    data=json_data,
                    file_name=f"troupeau_{date.today()}.json",
                    mime="application/json"
                )

# ============================================================================
# SECTION 12: PAGE PRODUCTION
# ============================================================================
def page_production():
    """Page de suivi de production"""
    st.markdown('<h2 class="section-header">🥛 SUIVI DE PRODUCTION</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📝 SAISIE", "📈 ANALYSE", "🏆 CLASSEMENT"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom, race FROM brebis WHERE sexe = 'F' AND statut = 'active'")
        femelles = cursor.fetchall()
        
        if femelles:
            with st.form("form_production"):
                brebis_sel = st.selectbox("Sélectionner une brebis", 
                                        [f"{f[1]} ({f[2]})" for f in femelles])
                
                col_prod1, col_prod2, col_prod3 = st.columns(3)
                
                with col_prod1:
                    quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5, 0.1)
                    cellules = st.number_input("Cellules (x1000)", 0, 1000, 200)
                
                with col_prod2:
                    mg = st.number_input("Matière grasse %", 0.0, 20.0, 7.2, 0.1)
                    lactose = st.number_input("Lactose %", 0.0, 10.0, 4.8, 0.1)
                
                with col_prod3:
                    proteine = st.number_input("Protéine %", 0.0, 20.0, 5.5, 0.1)
                    ph = st.number_input("pH", 6.0, 7.0, 6.7, 0.1)
                
                date_mesure = st.date_input("Date", value=date.today())
                notes = st.text_area("Notes")
                
                if st.form_submit_button("💾 Enregistrer", type="primary"):
                    try:
                        brebis_id = femelles[[f[1] for f in femelles].index(brebis_sel.split(' (')[0])][0]
                        
                        cursor.execute('''
                            INSERT INTO production_lait 
                            (brebis_id, date_mesure, quantite_litre, taux_matiere_grasse, 
                             taux_proteine, cellules_somatiques, lactose, ph, notes)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (brebis_id, date_mesure.isoformat(), quantite, mg, proteine, 
                              cellules*1000, lactose, ph, notes))
                        conn.commit()
                        st.success("✅ Production enregistrée!")
                    except:
                        st.error("Erreur lors de l'enregistrement")
        else:
            st.warning("Aucune brebis femelle active")
    
    with tab2:
        cursor.execute("""
            SELECT b.race,
                   strftime('%Y-%m', p.date_mesure) as mois,
                   AVG(p.quantite_litre) as lait_moyen,
                   AVG(p.taux_matiere_grasse) as mg_moyen,
                   AVG(p.taux_proteine) as proteine_moyenne,
                   COUNT(*) as mesures
            FROM production_lait p
            JOIN brebis b ON p.brebis_id = b.id
            WHERE p.date_mesure > date('now', '-6 months')
            GROUP BY b.race, mois
            HAVING mesures >= 2
            ORDER BY b.race, mois
        """)
        
        prod_data = cursor.fetchall()
        
        if prod_data:
            df_prod = pd.DataFrame(prod_data, 
                                 columns=['Race', 'Mois', 'Lait (L)', 'MG (%)', 'Protéine (%)', 'Mesures'])
            
            # Graphique évolution
            fig = px.line(df_prod, x='Mois', y='Lait (L)', color='Race',
                         title="Évolution production laitière par race",
                         markers=True)
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 🏆 CLASSEMENT DES PRODUCTRICES")
        
        cursor.execute("""
            SELECT b.nom, b.race, b.identifiant,
                   AVG(p.quantite_litre) as moyenne_lait,
                   AVG(p.taux_matiere_grasse) as moyenne_mg,
                   AVG(p.taux_proteine) as moyenne_proteine,
                   COUNT(*) as nb_mesures
            FROM production_lait p
            JOIN brebis b ON p.brebis_id = b.id
            WHERE p.date_mesure > date('now', '-90 days')
            GROUP BY b.id
            HAVING nb_mesures >= 3
            ORDER BY moyenne_lait DESC
            LIMIT 10
        """)
        
        top_prod = cursor.fetchall()
        
        if top_prod:
            df_top = pd.DataFrame(top_prod, 
                                columns=['Nom', 'Race', 'ID', 'Lait moyen (L)', 'MG (%)', 'Protéine (%)', 'Mesures'])
            
            st.dataframe(df_top.style.highlight_max(subset=['Lait moyen (L)'], color='lightgreen'))
            
            # Graphique top 5
            fig = px.bar(df_top.head(5), x='Nom', y='Lait moyen (L)',
                        color='Race',
                        title="Top 5 productrices - 90 derniers jours",
                        hover_data=['MG (%)', 'Protéine (%)', 'Mesures'])
            st.plotly_chart(fig, use_container_width=True)

# ============================================================================
# SECTION 13: PAGE CRITÈRES DE SÉLECTION
# ============================================================================
def page_criteres():
    """Page des critères de sélection morphologiques et phénotypiques"""
    st.markdown('<h2 class="section-header">🎯 CRITÈRES DE SÉLECTION - MAMMELLES</h2>', unsafe_allow_html=True)
    
    # Sélection de la brebis
    cursor = conn.cursor()
    cursor.execute("SELECT id, identifiant, nom, race, sexe FROM brebis WHERE sexe = 'F' ORDER BY race, identifiant")
    brebis_femelles = cursor.fetchall()
    
    if not brebis_femelles:
        st.warning("Aucune brebis femelle dans la base de données")
        return
    
    brebis_option = st.selectbox(
        "SÉLECTIONNEZ UNE BREBIS FAMELLE:",
        [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_femelles],
        key="criteres_selection"
    )
    
    if not brebis_option:
        return
    
    try:
        brebis_id = int(brebis_option.split('(')[1].split(')')[0].split('-')[-1])
    except:
        st.error("Erreur dans le format de l'identifiant")
        return
    
    # Récupérer les données de la brebis
    cursor.execute("SELECT * FROM brebis WHERE id = ?", (brebis_id,))
    brebis_info = cursor.fetchone()
    
    if not brebis_info:
        st.error("Brebis non trouvée")
        return
    
    columns = [desc[0] for desc in cursor.description]
    brebis_dict = dict(zip(columns, brebis_info))
    
    tab1, tab2, tab3, tab4 = st.tabs(["📏 MAMMELLES", "🏋️ MORPHOLOGIE", "🧬 PHÉNOTYPE", "📊 SCORING"])
    
    with tab1:
        st.markdown("### 📏 CRITÈRES MAMMAIRES - PRODUCTION LAITIÈRE")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class='mammelle-card'>
                <h4>🎯 PRINCIPAUX CRITÈRES</h4>
                <p><strong>1. Volume mammaire</strong> (1-5)</p>
                <p><strong>2. Symétrie</strong> (1-5)</p>
                <p><strong>3. Insertion des trayons</strong> (1-5)</p>
                <p><strong>4. Longueur des trayons</strong> (cm)</p>
                <p><strong>5. Orientation</strong></p>
            </div>
            """, unsafe_allow_html=True)
            
            # Afficher les critères existants
            st.markdown("#### 📋 CRITÈRES ACTUELS")
            if brebis_dict.get('volume_mammaire'):
                col_crit1, col_crit2, col_crit3 = st.columns(3)
                with col_crit1:
                    st.metric("Volume", f"{brebis_dict['volume_mammaire']}/5")
                with col_crit2:
                    st.metric("Symétrie", f"{brebis_dict['symetrie_mammaire']}/5")
                with col_crit3:
                    st.metric("Insertion", f"{brebis_dict['insertion_trayons']}/5")
                
                if brebis_dict.get('longueur_trayons_cm'):
                    st.metric("Longueur trayons", f"{brebis_dict['longueur_trayons_cm']:.1f} cm")
            
            # Formulaire d'évaluation des mamelles
            with st.form("evaluation_mamelles"):
                st.markdown("#### 📝 ÉVALUATION MAMMAIRE")
                
                volume = st.slider("Volume mammaire (1-5)", 1, 5, 
                                  value=brebis_dict.get('volume_mammaire', 3) or 3,
                                  help="1: Très petit, 5: Très développé")
                symetrie = st.slider("Symétrie (1-5)", 1, 5,
                                    value=brebis_dict.get('symetrie_mammaire', 3) or 3,
                                    help="1: Asymétrique, 5: Parfaitement symétrique")
                insertion = st.slider("Insertion des trayons (1-5)", 1, 5,
                                     value=brebis_dict.get('insertion_trayons', 3) or 3,
                                     help="1: Très écartés, 5: Bien insérés")
                longueur_trayons = st.slider("Longueur des trayons (cm)", 2.0, 8.0, 
                                          value=brebis_dict.get('longueur_trayons_cm', 4.5) or 4.5, step=0.1
                orientation = st.selectbox("Orientation des trayons",
                                         ['parallele', 'leger_divergent', 'divergent'],
                                         index=['parallele', 'leger_divergent', 'divergent'].index(
                                             brebis_dict.get('orientation_trayons', 'parallele') or 'parallele'
                                         ))
                
                if st.form_submit_button("💾 Enregistrer évaluation", type="primary"):
                    cursor.execute('''
                        UPDATE brebis 
                        SET volume_mammaire = ?, symetrie_mammaire = ?, insertion_trayons = ?,
                            longueur_trayons_cm = ?, orientation_trayons = ?
                        WHERE id = ?
                    ''', (volume, symetrie, insertion, longueur_trayons, orientation, brebis_id))
                    conn.commit()
                    
                    score_total = (volume + symetrie + insertion) / 3
                    st.success(f"✅ Évaluation enregistrée! Score mammelle: {score_total:.1f}/5")
        
        with col2:
            # Standards par race
            race = brebis_dict.get('race', 'INCONNU')
            if race in STANDARDS_RACES and 'criteres_mammaires' in STANDARDS_RACES[race]:
                criteres = STANDARDS_RACES[race]['criteres_mammaires']
                st.markdown(f"""
                <div class='race-card'>
                    <h4>🏷️ STANDARDS {race}</h4>
                    <p><strong>Volume:</strong> {criteres['volume']}</p>
                    <p><strong>Trayons:</strong> {criteres['trayons']}</p>
                    <p><strong>Symétrie:</strong> {criteres['symetrie']}</p>
                    <p><strong>Aptitude laitière:</strong> {criteres['aptitude_laitiere']}</p>
                </div>
                """, unsafe_allow_html=True)
            
            # Classification mammaire
            st.markdown("""
            <div style='text-align: center; padding: 20px; background: #f8f9fa; border-radius: 10px;'>
                <h4>📊 CLASSIFICATION MAMMAIRE</h4>
                <p><strong>Type A (4-5):</strong> Excellent pour la production</p>
                <p><strong>Type B (3-4):</strong> Bon pour la production</p>
                <p><strong>Type C (2-3):</strong> Moyen, à surveiller</p>
                <p><strong>Type D (1-2):</strong> À améliorer ou réformer</p>
            </div>
            """, unsafe_allow_html=True)
    
    with tab2:
        st.markdown("### 🏋️ CRITÈRES MORPHOLOGIQUES GÉNÉRAUX")
        
        col_morph1, col_morph2 = st.columns(2)
        
        with col_morph1:
            # Mensurations actuelles
            st.markdown("#### 📏 MENSURATIONS ACTUELLES")
            if brebis_dict.get('longueur_corps_cm'):
                col_meas1, col_meas2 = st.columns(2)
                with col_meas1:
                    st.metric("Longueur", f"{brebis_dict['longueur_corps_cm']:.1f} cm")
                    st.metric("Hauteur", f"{brebis_dict['hauteur_garrot_cm']:.1f} cm")
                with col_meas2:
                    st.metric("Largeur bassin", f"{brebis_dict['largeur_bassin_cm']:.1f} cm")
                    st.metric("Tour poitrine", f"{brebis_dict['tour_poitrine_cm']:.1f} cm")
            
            # Standards de race
            if race in STANDARDS_RACES:
                standards = STANDARDS_RACES[race]['mensurations']
                st.markdown("#### 🎯 STANDARDS DE RACE")
                st.write(f"**Longueur:** {standards['longueur_cm'][0]}-{standards['longueur_cm'][1]} cm")
                st.write(f"**Hauteur:** {standards['hauteur_cm'][0]}-{standards['hauteur_cm'][1]} cm")
                st.write(f"**Largeur bassin:** {standards['largeur_bassin_cm'][0]}-{standards['largeur_bassin_cm'][1]} cm")
        
        with col_morph2:
            # Caractères qualitatifs
            st.markdown("#### 🎨 CARACTÈRES QUALITATIFS")
            
            caract_data = {
                'Caractère': ['Couleur', 'Intensité couleur', 'Cornes', 'Type laine', 'Qualité laine', 'Tempérament'],
                'Valeur': [
                    brebis_dict.get('couleur_robe', 'N/A'),
                    f"{brebis_dict.get('intensite_couleur', 0)}/10",
                    'Oui' if brebis_dict.get('cornes') else 'Non',
                    brebis_dict.get('type_laine', 'N/A'),
                    f"{brebis_dict.get('qualite_laine', 0)}/10",
                    brebis_dict.get('temperement', 'N/A')
                ]
            }
            
            df_caract = pd.DataFrame(caract_data)
            st.dataframe(df_caract, use_container_width=True, hide_index=True)
    
    with tab3:
        st.markdown("### 🧬 CRITÈRES PHÉNOTYPIQUES")
        
        # Phénotypes liés à la production
        st.markdown("#### 🥛 PHÉNOTYPES LAITIERS")
        
        phenotypes = pd.DataFrame({
            'Caractère': ['Développement mammaire', 'Tempérament', 'Appétit', 
                         'Rapidité de traite', 'Résistance mammite'],
            'Héritabilité': [0.35, 0.15, 0.20, 0.30, 0.10],
            'Impact production': ['Élevé', 'Faible', 'Moyen', 'Élevé', 'Élevé']
        })
        
        st.dataframe(phenotypes.style.background_gradient(subset=['Héritabilité']))
    
    with tab4:
        st.markdown("### 📊 SCORING INTÉGRÉ")
        
        # Formulaire complet de scoring
        with st.form("scoring_complet"):
            st.markdown("#### 🎯 ÉVALUATION COMPLÈTE")
            
            # Récupérer les valeurs actuelles
            vol_mam = brebis_dict.get('volume_mammaire', 3) or 3
            sym_mam = brebis_dict.get('symetrie_mammaire', 3) or 3
            ins_mam = brebis_dict.get('insertion_trayons', 3) or 3
            
            col_score1, col_score2, col_score3 = st.columns(3)
            
            with col_score1:
                st.markdown("**📏 MORPHOLOGIE**")
                conformation = st.slider("Conformation générale (1-10)", 1, 10, 
                                       value=int(brebis_dict.get('score_conformation', 7) or 7))
                developpement = st.slider("Développement musculaire (1-10)", 1, 10, 6)
            
            with col_score2:
                st.markdown("**🥛 MAMMELLES**")
                volume_m = st.slider("Volume mammaire (1-10)", 1, 10, vol_mam * 2)
                symetrie_m = st.slider("Symétrie mammaire (1-10)", 1, 10, sym_mam * 2)
                insertion_m = st.slider("Insertion trayons (1-10)", 1, 10, ins_mam * 2)
            
            with col_score3:
                st.markdown("**🧬 GÉNÉTIQUE**")
                valeur_gen = st.slider("Valeur génétique estimée (1-10)", 1, 10, 7)
                diversite = st.slider("Diversité génétique (1-10)", 1, 10, 8)
            
            notes = st.text_area("Observations complémentaires", 
                                value=brebis_dict.get('notes', ''))
            
            if st.form_submit_button("🎯 Calculer score final", type="primary"):
                # Calcul des scores
                score_morph = (conformation + developpement) / 2
                score_mamelle = (volume_m + symetrie_m + insertion_m) / 3
                score_gen = (valeur_gen + diversite) / 2
                score_final = (score_morph * 0.3 + score_mamelle * 0.4 + score_gen * 0.3)
                
                # Mettre à jour la base
                cursor.execute('''
                    UPDATE brebis 
                    SET score_conformation = ?, aptitudes = ?
                    WHERE id = ?
                ''', (score_final, notes, brebis_id))
                conn.commit()
                
                # Affichage résultats
                col_res1, col_res2, col_res3 = st.columns(3)
                
                with col_res1:
                    st.metric("Score morphologie", f"{score_morph:.1f}/10")
                
                with col_res2:
                    st.metric("Score mamelles", f"{score_mamelle:.1f}/10")
                
                with col_res3:
                    st.metric("Score génétique", f"{score_gen:.1f}/10")
                
                # Score final avec interprétation
                st.markdown("---")
                st.markdown(f"### 🏆 SCORE FINAL: **{score_final:.1f}/10**")
                
                if score_final >= 8:
                    st.success("🎖️ **EXCELLENT** - Animal d'élite pour la reproduction")
                elif score_final >= 6:
                    st.info("✅ **BON** - Animal de production satisfaisant")
                elif score_final >= 4:
                    st.warning("⚠️ **MOYEN** - À surveiller ou améliorer")
                else:
                    st.error("❌ **FAIBLE** - À réformer ou surveiller étroitement")

# ============================================================================
# SECTION 14: PAGE STATISTIQUES (RSTATS)
# ============================================================================
def page_stats():
    """Page d'analyse statistique avancée"""
    st.markdown('<h2 class="section-header">📊 ANALYSE STATISTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    cursor = conn.cursor()
    cursor.execute("""
        SELECT race, sexe, age_mois, poids, score_condition, 
               longueur_corps_cm, hauteur_garrot_cm, largeur_bassin_cm,
               tour_poitrine_cm, intensite_couleur
        FROM brebis
        WHERE statut = 'active'
    """)
    
    data = cursor.fetchall()
    df = pd.DataFrame(data, columns=['Race', 'Sexe', 'Âge', 'Poids', 'Score', 'Longueur', 
                                     'Hauteur', 'Largeur', 'Poitrine', 'Intensité'])
    
    tab1, tab2 = st.tabs(["📈 DESCRIPTIVE", "📊 CORRÉLATIONS"])
    
    with tab1:
        st.markdown("### 📈 ANALYSE DESCRIPTIVE PAR RACE")
        
        races = df['Race'].unique()
        for race in races:
            with st.expander(f"{STANDARDS_RACES.get(race, {}).get('nom_complet', race)}"):
                race_df = df[df['Race'] == race]
                
                col_stat1, col_stat2, col_stat3 = st.columns(3)
                
                with col_stat1:
                    st.metric("Nombre", len(race_df))
                    st.metric("Poids moyen", f"{race_df['Poids'].mean():.1f} kg")
                
                with col_stat2:
                    st.metric("Âge moyen", f"{race_df['Âge'].mean():.0f} mois")
                    st.metric("Score moyen", f"{race_df['Score'].mean():.1f}/5")
                
                with col_stat3:
                    st.metric("Longueur moyenne", f"{race_df['Longueur'].mean():.1f} cm")
                    st.metric("Hauteur moyenne", f"{race_df['Hauteur'].mean():.1f} cm")
    
    with tab2:
        st.markdown("### 📊 MATRICE DE CORRÉLATION")
        
        numeric_vars = ['Poids', 'Longueur', 'Hauteur', 'Largeur', 'Âge', 'Score']
        numeric_df = df[numeric_vars].dropna()
        
        if not numeric_df.empty:
            corr_matrix = numeric_df.corr()
            
            fig = px.imshow(corr_matrix,
                           title="Matrice de corrélation - Variables morphologiques",
                           color_continuous_scale='RdBu',
                           zmin=-1, zmax=1,
                           text_auto=True)
            st.plotly_chart(fig, use_container_width=True)

# ============================================================================
# SECTION 15: PAGE GÉNÉTIQUE
# ============================================================================
def page_genetique():
    """Page d'analyse génétique avancée"""
    st.markdown('<h2 class="section-header">🧬 ANALYSE GÉNÉTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🧪 GÉNOTYPAGE", "🌳 DIVERSITÉ"])
    
    with tab1:
        st.markdown("### 🧪 GÉNOTYPAGE SNP")
        
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY race")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_select = st.selectbox("Sélectionner une brebis", 
                                        [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list],
                                        key="geno_selection")
            
            if brebis_select and st.button("🧬 Générer génotype", type="primary", key="geno_btn"):
                try:
                    brebis_id = int(brebis_select.split('(')[1].split(')')[0].split('-')[-1])
                    race = brebis_select.split('- ')[1]
                    
                    genotypes = ModuleGenetique.generer_genotype(brebis_id, race)
                    
                    # Insérer dans la base
                    for genotype in genotypes:
                        cursor.execute('''
                            INSERT OR REPLACE INTO genotypage 
                            (brebis_id, marqueur, chromosome, position, allele1, allele2, 
                             genotype, frequence_allelique, effet_additif, effet_dominant, 
                             r2, p_value, gene_associe, trait_associe, date_analyse)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', genotype)
                    
                    conn.commit()
                    st.success(f"✅ Génotype généré pour {brebis_select}")
                    
                except Exception as e:
                    st.error(f"Erreur: {str(e)}")
    
    with tab2:
        st.markdown("### 🌳 DIVERSITÉ GÉNÉTIQUE")
        
        cursor.execute("""
            SELECT g.brebis_id, b.race, b.identifiant,
                   g.marqueur, g.allele1, g.allele2, g.genotype,
                   g.frequence_allelique, g.trait_associe
            FROM genotypage g
            JOIN brebis b ON g.brebis_id = b.id
            LIMIT 50
        """)
        
        geno_data = cursor.fetchall()
        
        if geno_data:
            try:
                diversite = ModuleGenetique.calculer_diversite_genetique(geno_data)
                
                if diversite:
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        st.metric("Hétérozygotie observée", f"{diversite['heterozygosite_observee']:.4f}")
                    with col2:
                        st.metric("Hétérozygotie attendue", f"{diversite['heterozygosite_attendue']:.4f}")
                    with col3:
                        st.metric("Fis", f"{diversite['fis']:.4f}")
                    with col4:
                        st.metric("SNPs analysés", diversite['nombre_snps'])
                else:
                    st.info("Impossible de calculer la diversité génétique")
                    
            except Exception as e:
                st.error(f"Erreur dans le calcul: {str(e)}")
        else:
            st.info("Aucune donnée de génotypage disponible. Générez d'abord des génotypes.")

# ============================================================================
# SECTION 16: BARRE LATÉRALE
# ============================================================================
with st.sidebar:
    st.markdown("""
    <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #1a237e 0%, #283593 100%); 
                color: white; border-radius: 10px; margin-bottom: 20px;'>
        <h2>🐑 RACES ALGÉRIENNES</h2>
        <p>Système de gestion scientifique</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("### 📍 NAVIGATION")
    
    page = st.radio(
        "MENU PRINCIPAL",
        ["🏠 ACCUEIL", 
         "📐 SCANNER 3D", 
         "📊 GESTION", 
         "🥛 PRODUCTION",
         "🎯 CRITÈRES",
         "📊 RSTATS",
         "🧬 GÉNÉTIQUE"]
    )
    
    st.markdown("---")
    
    # Statistiques rapides
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE statut = 'active'")
    actives = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F' AND statut = 'active'")
    femelles = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genotypage")
    genotypages = cursor.fetchone()[0]
    
    st.markdown("### 📊 EN DIRECT")
    st.metric("🐑 Actives", actives)
    st.metric("♀️ Femelles", femelles)
    st.metric("🧬 Génotypages", genotypages)
    
    st.markdown("---")
    
    # Standards des races
    st.markdown("### 🏷️ STANDARDS RACES")
    
    race_info = st.selectbox("Info race", list(STANDARDS_RACES.keys()))
    
    if race_info in STANDARDS_RACES:
        info = STANDARDS_RACES[race_info]
        st.markdown(f"""
        <div class='race-card' style='padding: 10px;'>
            <h5>{info['nom_complet']}</h5>
            <p><small>Poids ♀️: {info['poids_adulte']['femelle'][0]}-{info['poids_adulte']['femelle'][1]} kg</small></p>
            <p><small>Poids ♂️: {info['poids_adulte']['male'][0]}-{info['poids_adulte']['male'][1]} kg</small></p>
            <p><small>Production: {info['production_lait'][0]}-{info['production_lait'][1]} L/j</small></p>
        </div>
        """, unsafe_allow_html=True)

# ============================================================================
# SECTION 17: NAVIGATION PRINCIPALE
# ============================================================================
if page == "🏠 ACCUEIL":
    page_accueil()
elif page == "📐 SCANNER 3D":
    page_scanner_3d()
elif page == "📊 GESTION":
    page_gestion()
elif page == "🥛 PRODUCTION":
    page_production()
elif page == "🎯 CRITÈRES":
    page_criteres()
elif page == "📊 RSTATS":
    page_stats()
elif page == "🧬 GÉNÉTIQUE":
    page_genetique()

# ============================================================================
# SECTION 18: PIED DE PAGE
# ============================================================================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version 4.0</p>
    <p>📐 Scanner 3D • 🎯 Critères de sélection • 🧬 Génétique • 📊 Statistiques</p>
    <p>© 2024 - Système de gestion scientifique des races ovines algériennes</p>
</div>
""", unsafe_allow_html=True)
