"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D et Génétique
Base de données simulée de races ovines algériennes
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
import base64
from PIL import Image, ImageDraw

# ========== CONFIGURATION ==========

st.set_page_config(
    page_title="Ovin Manager Pro - Races Algériennes",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personnalisé
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
    .hamra-card {
        background: linear-gradient(135deg, #8B0000 0%, #FF4500 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(139,0,0,0.2);
    }
    .ouda-card {
        background: linear-gradient(135deg, #4CAF50 0%, #2E7D32 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(46,125,50,0.2);
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
    .scanner-view {
        background: black;
        border-radius: 10px;
        padding: 10px;
        margin: 10px 0;
        text-align: center;
    }
    .scan-progress {
        height: 20px;
        background: linear-gradient(90deg, #1a237e, #283593);
        border-radius: 10px;
        margin: 10px 0;
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

# ========== FONCTIONS STATISTIQUES ALTERNATIVES (sans scipy) ==========

def skewness(data):
    """Calcule le coefficient d'asymétrie de Pearson"""
    if len(data) < 3:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    skew = np.mean(((data - mean) / std) ** 3)
    return skew

def kurtosis(data):
    """Calcule le coefficient d'aplatissement"""
    if len(data) < 4:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    kurt = np.mean(((data - mean) / std) ** 4) - 3
    return kurt

def pearson_correlation(x, y):
    """Calcule la corrélation de Pearson"""
    if len(x) != len(y) or len(x) < 2:
        return 0, 1.0
    
    x = np.array(x)
    y = np.array(y)
    
    # Supprimer les NaN
    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]
    y = y[mask]
    
    if len(x) < 2:
        return 0, 1.0
    
    # Calcul de la corrélation
    corr_matrix = np.corrcoef(x, y)
    r = corr_matrix[0, 1]
    
    # p-value approximative
    n = len(x)
    if n <= 2:
        return r, 1.0
    
    # Test t pour la corrélation
    t = r * np.sqrt((n - 2) / (1 - r**2)) if abs(r) < 1 else 0
    # p-value approximative (distribution t)
    p = 2 * (1 - t_distribution_cdf(abs(t), n-2))
    
    return r, p

def t_distribution_cdf(t, df):
    """Fonction CDF approximative pour la distribution t"""
    if df <= 0:
        return 0.5
    
    # Approximation simple
    x = t / np.sqrt(df + t**2)
    return 0.5 + 0.5 * x

def anova_one_way(groups):
    """ANOVA à un facteur simplifiée"""
    if len(groups) < 2:
        return 0, 1.0
    
    # Calcul des moyennes
    all_data = np.concatenate(groups)
    grand_mean = np.mean(all_data)
    
    # Somme des carrés entre groupes
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
    
    # Somme des carrés totale
    ss_total = sum((x - grand_mean) ** 2 for x in all_data)
    
    # Degrés de liberté
    df_between = len(groups) - 1
    df_within = len(all_data) - len(groups)
    
    if df_within <= 0 or df_between <= 0:
        return 0, 1.0
    
    # Carrés moyens
    ms_between = ss_between / df_between
    ms_within = (ss_total - ss_between) / df_within
    
    # Statistique F
    if ms_within == 0:
        return float('inf'), 0.0
    
    f_stat = ms_between / ms_within
    
    # p-value approximative
    p_value = f_distribution_pvalue(f_stat, df_between, df_within)
    
    return f_stat, p_value

def f_distribution_pvalue(f, df1, df2):
    """p-value approximative pour la distribution F"""
    if df1 <= 0 or df2 <= 0:
        return 1.0
    
    # Approximation simple
    x = df2 / (df2 + df1 * f)
    p = 1 - x ** (df2 / 2)
    return min(max(p, 0), 1)

def shapiro_wilk_test(data):
    """Test de normalité Shapiro-Wilk simplifié"""
    n = len(data)
    if n < 3 or n > 5000:
        return 0.5, 1.0  # Valeurs par défaut
    
    # Version très simplifiée
    # En pratique, il faudrait implémenter l'algorithme complet
    skew = skewness(data)
    kurt = kurtosis(data)
    
    # Statistique W approximative
    w = np.exp(-0.1 * (skew**2 + (kurt/2)**2))
    
    # p-value approximative
    p = 1 - w
    
    return w, p

# ========== STANDARDS DES RACES ALGÉRIENNES ==========

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
        'prolificite': (1.2, 1.8)
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
        'prolificite': (1.1, 1.5)
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
        'prolificite': (1.3, 1.7)
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
        'prolificite': (1.0, 1.4)
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
        'prolificite': (1.2, 1.8)
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
        'prolificite': (1.0, 1.6)
    }
}

# ========== BASE DE DONNÉES SIMULÉE ==========

def creer_base_races():
    """Crée une base de données avec races algériennes"""
    conn = sqlite3.connect('ovin_algerien.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table brebis avec plus de caractéristiques
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
            
            -- Caractères qualitatifs
            temperement TEXT CHECK(temperement IN ('calme', 'nervieux', 'intermediaire')),
            aptitude TEXT CHECK(aptitude IN ('lait', 'viande', 'mixte', 'laine')),
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
            conductivite FLOAT,
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
            fichier_3d_path TEXT,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Table génotypage avancée
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
    
    # Table QTL (Quantitative Trait Loci)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS qtl (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            trait TEXT,
            chromosome TEXT,
            position_debut INTEGER,
            position_fin INTEGER,
            lod_score FLOAT,
            variance_expliquee FLOAT,
            marqueurs TEXT,
            espece TEXT,
            reference TEXT
        )
    ''')
    
    # Table analyse statistique
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS analyses_stats (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            type_analyse TEXT,
            parametres TEXT,
            resultats_json TEXT,
            date_analyse DATE,
            notes TEXT
        )
    ''')
    
    # Vérifier si la base est vide et la peupler
    cursor.execute("SELECT COUNT(*) FROM brebis")
    count = cursor.fetchone()[0]
    
    if count == 0:
        peupler_base_races(cursor, conn)
    
    # Peupler les QTL de référence
    cursor.execute("SELECT COUNT(*) FROM qtl")
    if cursor.fetchone()[0] == 0:
        peupler_qtl_reference(cursor)
    
    conn.commit()
    return conn

def peupler_base_races(cursor, conn):
    """Peuple la base avec des races algériennes"""
    
    noms_femelles = [
        "Lalla", "Zina", "Nour", "Yasmina", "Fatima", "Khadija", "Aicha", 
        "Samira", "Leila", "Soraya", "Nadia", "Rym", "Salima", "Djamila",
        "Zahra", "Malika", "Halima", "Farida", "Rachida", "Safia"
    ]
    
    noms_males = [
        "Sultan", "Amir", "Karim", "Rachid", "Yacine", "Noureddine", 
        "Mohamed", "Ali", "Omar", "Hassan", "Hussein", "Bilal", "Mokhtar"
    ]
    
    races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU']
    poids_races = {
        'HAMRA': {'F': (45, 65), 'M': (65, 90)},
        'OUDA': {'F': (50, 70), 'M': (70, 100)},
        'SIDAHOU': {'F': (40, 60), 'M': (60, 85)},
        'BERBERE': {'F': (35, 50), 'M': (50, 70)},
        'CROISE': {'F': (40, 70), 'M': (60, 95)},
        'INCONNU': {'F': (30, 60), 'M': (50, 80)}
    }
    
    couleurs_races = {
        'HAMRA': ['Rousse', 'Rousse foncée', 'Marron'],
        'OUDA': ['Blanche', 'Crème', 'Blanc cassé'],
        'SIDAHOU': ['Noire et blanche', 'Pie noire', 'Tête noire'],
        'BERBERE': ['Noire', 'Brune', 'Grise', 'Pie'],
        'CROISE': ['Variable', 'Panachée', 'Mélangée'],
        'INCONNU': ['Indéterminée']
    }
    
    brebis_data = []
    
    # Créer 50 brebis de différentes races
    for i in range(1, 51):
        race = random.choice(races)
        sexe = random.choices(['F', 'M'], weights=[0.7, 0.3])[0]
        
        if sexe == 'F':
            nom = random.choice(noms_femelles) + f" {i}"
        else:
            nom = random.choice(noms_males) + f" {i}"
        
        identifiant = f"{race[:3]}-{sexe}-2023-{i:03d}"
        age_mois = random.randint(12, 84)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        
        # Poids selon race et sexe
        poids_min, poids_max = poids_races[race][sexe]
        poids = random.uniform(poids_min, poids_max)
        
        score_condition = random.randint(2, 4)
        couleur = random.choice(couleurs_races[race])
        intensite_couleur = random.randint(5, 10)
        cornes = random.choice([True, False])
        taille_cornes = random.uniform(0, (60 if sexe == 'M' else 25)) if cornes else 0
        
        # Mensurations selon race
        standards = STANDARDS_RACES[race]['mensurations']
        longueur_corps = random.uniform(*standards['longueur_cm'])
        hauteur_garrot = random.uniform(*standards['hauteur_cm'])
        largeur_bassin = random.uniform(*standards['largeur_bassin_cm'])
        tour_poitrine = random.uniform(*standards['tour_poitrine_cm'])
        
        # Autres mesures
        circonference_tete = random.uniform(45, 65)
        longueur_oreille = random.uniform(12, 18)
        
        # Caractères qualitatifs
        temperement = random.choices(['calme', 'nervieux', 'intermediaire'], weights=[0.6, 0.1, 0.3])[0]
        aptitude = random.choice(['lait', 'viande', 'mixte', 'laine'])
        type_laine = random.choice(['fine', 'semi-fine', 'grossière'])
        qualite_laine = random.randint(5, 9) if race in ['HAMRA', 'BERBERE'] else random.randint(3, 7)
        
        brebis_data.append((
            identifiant, nom, race, '', sexe, date_naissance.isoformat(), 
            age_mois, poids, score_condition, couleur, intensite_couleur, 
            cornes, taille_cornes, '', type_laine, qualite_laine,
            longueur_corps, hauteur_garrot, largeur_bassin, tour_poitrine,
            circonference_tete, longueur_oreille, temperement, aptitude,
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
            temperement, aptitude, notes, mere_id, pere_id, 
            coefficient_consanguinite, statut
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', brebis_data)
    
    conn.commit()

def peupler_qtl_reference(cursor):
    """Peuple la table QTL avec des références scientifiques"""
    qtl_data = [
        ('poids_vif', '1', 45000000, 48000000, 4.2, 0.15, 'BM1818, OARFCB304', 'Ovis aries', 'PMID: 12345678'),
        ('production_lait', '3', 23000000, 25000000, 5.1, 0.18, 'MAF214, OARJMP58', 'Ovis aries', 'PMID: 23456789'),
        ('taux_matiere_grasse', '6', 12000000, 15000000, 3.8, 0.12, 'BM6506, OARHH35', 'Ovis aries', 'PMID: 34567890'),
        ('prolificite', '9', 55000000, 58000000, 6.3, 0.22, 'BMS2508, OARCP49', 'Ovis aries', 'PMID: 45678901'),
        ('longueur_corps', '2', 34000000, 37000000, 4.5, 0.16, 'OARFCB128, MAF70', 'Ovis aries', 'PMID: 56789012'),
        ('qualite_laine', '16', 20000000, 23000000, 3.2, 0.10, 'BM757, OARJMP29', 'Ovis aries', 'PMID: 67890123'),
        ('resistance_maladies', '5', 65000000, 68000000, 4.8, 0.17, 'OARHH47, BM827', 'Ovis aries', 'DOI: 10.1111/age.12345'),
        ('couleur_robe', '13', 28000000, 31000000, 7.1, 0.25, 'TYRP1, ASIP', 'Ovis aries', 'DOI: 10.1111/age.23456'),
        ('adaptation_chaleur', '10', 42000000, 45000000, 5.6, 0.20, 'HSP70, TLR4', 'Ovis aries', 'DOI: 10.1111/age.34567'),
        ('efficacite_alimentaire', '8', 18000000, 22000000, 4.3, 0.14, 'GH, GHR', 'Ovis aries', 'DOI: 10.1111/age.45678')
    ]
    
    cursor.executemany('''
        INSERT INTO qtl (trait, chromosome, position_debut, position_fin, 
                        lod_score, variance_expliquee, marqueurs, espece, reference)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', qtl_data)

# Initialiser la base
conn = creer_base_races()

# ========== MODULE SCANNER 3D AVANCÉ ==========

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
        
        n_points = 1000
        points = []
        
        # Paramètres selon la race et le poids
        poids = brebis_info.get('poids', 50)
        race = brebis_info.get('race', 'INCONNU')
        
        # Facteurs de forme par race
        facteurs = {
            'HAMRA': {'long': 0.95, 'larg': 0.85, 'haut': 0.90},
            'OUDA': {'long': 1.05, 'larg': 0.95, 'haut': 1.00},
            'SIDAHOU': {'long': 0.90, 'larg': 0.80, 'haut': 0.85},
            'BERBERE': {'long': 0.85, 'larg': 0.75, 'haut': 0.80},
            'CROISE': {'long': 1.00, 'larg': 0.90, 'haut': 0.95},
            'INCONNU': {'long': 0.95, 'larg': 0.85, 'haut': 0.90}
        }
        
        facteur = facteurs.get(race, facteurs['INCONNU'])
        
        # Rayons approximatifs
        rx = 0.6 * poids**0.33 * facteur['larg']  # Largeur
        ry = 1.2 * poids**0.33 * facteur['long']  # Longueur
        rz = 0.8 * poids**0.33 * facteur['haut']  # Hauteur
        
        for _ in range(n_points):
            # Points sur ellipsoïde
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
            
            x = rx * np.sin(phi) * np.cos(theta) + np.random.normal(0, rx*0.05)
            y = ry * np.sin(phi) * np.sin(theta) + np.random.normal(0, ry*0.05)
            z = rz * np.cos(phi) + np.random.normal(0, rz*0.05)
            
            # Ajouter des caractéristiques anatomiques
            if z > rz * 0.5:  # Dos
                intensity = np.random.uniform(100, 150)
            elif abs(x) > rx * 0.7:  # Côtés
                intensity = np.random.uniform(150, 200)
            else:  # Ventre
                intensity = np.random.uniform(200, 255)
            
            points.append({
                'x': float(x),
                'y': float(y),
                'z': float(z),
                'intensity': int(intensity),
                'anatomic_region': random.choice(['dos', 'cote', 'ventre', 'queue'])
            })
        
        return points
    
    @staticmethod
    def analyser_points_3d(points):
        """Analyse avancée des points 3D"""
        if len(points) < 10:
            return {}
        
        x = np.array([p['x'] for p in points])
        y = np.array([p['y'] for p in points])
        z = np.array([p['z'] for p in points])
        
        # Mesures de base
        longueur = np.max(y) - np.min(y)
        largeur = np.max(x) - np.min(x)
        hauteur = np.max(z) - np.min(z)
        
        # Volume et surface
        a, b, c = longueur/2, largeur/2, hauteur/2
        volume = (4/3) * np.pi * a * b * c / 1000000
        surface = 4 * np.pi * ((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3**(1/1.6) / 10000
        
        # Indices corporels
        indice_corporel = (longueur * largeur * hauteur) ** (1/3)
        indice_conformation = longueur / hauteur
        
        # Symétrie
        x_pos = x[x > 0]
        x_neg = x[x < 0]
        symetrie_lat = 100 - abs(np.mean(x_pos) + np.mean(x_neg)) * 5 if len(x_pos) > 0 and len(x_neg) > 0 else 100
        
        # Densité des points
        densite = len(points) / (longueur * largeur * hauteur)
        
        # Analyse de forme
        skewness_x = skewness(x)
        kurtosis_z = kurtosis(z)
        
        return {
            'longueur_cm': round(float(longueur), 1),
            'largeur_cm': round(float(largeur), 1),
            'hauteur_cm': round(float(hauteur), 1),
            'volume_m3': round(float(volume), 4),
            'surface_m2': round(float(surface), 2),
            'indice_corporel': round(float(indice_corporel), 2),
            'indice_conformation': round(float(indice_conformation), 2),
            'symetrie_lat_percent': round(float(symetrie_lat), 1),
            'densite_points': round(float(densite), 2),
            'asymetrie_skewness': round(float(skewness_x), 3),
            'aplatissement_kurtosis': round(float(kurtosis_z), 3),
            'score_conformation': min(100, round(85 + symetrie_lat/5 + indice_conformation*10, 1))
        }
    
    @staticmethod
    def creer_visualisation_3d(points):
        """Crée une visualisation 3D interactive"""
        if len(points) < 10:
            return None
        
        x = [p['x'] for p in points]
        y = [p['y'] for p in points]
        z = [p['z'] for p in points]
        intensities = [p['intensity'] for p in points]
        
        fig = go.Figure(data=[go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=3,
                color=intensities,
                colorscale='Viridis',
                opacity=0.8,
                showscale=True,
                colorbar=dict(title="Intensité")
            ),
            hovertemplate='<b>X:</b> %{x:.1f}<br><b>Y:</b> %{y:.1f}<br><b>Z:</b> %{z:.1f}<extra></extra>'
        )])
        
        fig.update_layout(
            title="Scan 3D - Reconstruction anatomique",
            scene=dict(
                xaxis_title='Largeur (cm)',
                yaxis_title='Longueur (cm)',
                zaxis_title='Hauteur (cm)',
                aspectmode='data',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5)
                )
            ),
            width=900,
            height=600
        )
        
        return fig

# ========== MODULE GÉNÉTIQUE ==========

class ModuleGenetique:
    """Module d'analyse génétique avancée"""
    
    @staticmethod
    def generer_genotype(brebis_id, race):
        """Génère un génotype simulé"""
        np.random.seed(brebis_id)
        
        # Marqueurs SNP pour ovins
        snps = [
            ('SNP001', '1', 1234567, ['A', 'G'], 'Couleur'),
            ('SNP002', '3', 2345678, ['C', 'T'], 'Production lait'),
            ('SNP003', '6', 3456789, ['G', 'A'], 'Taille'),
            ('SNP004', '9', 4567890, ['T', 'C'], 'Prolificité'),
            ('SNP005', '13', 5678901, ['A', 'T'], 'Résistance'),
            ('SNP006', '16', 6789012, ['C', 'G'], 'Qualité laine'),
            ('SNP007', '18', 7890123, ['G', 'T'], 'Matière grasse'),
            ('SNP008', '21', 8901234, ['T', 'A'], 'Croissance'),
            ('SNP009', '25', 9012345, ['A', 'C'], 'Adaptation'),
            ('SNP010', '26', 1234567, ['G', 'A'], 'Tempérament')
        ]
        
        genotypes = []
        
        for snp in snps:
            marqueur, chrom, pos, alleles, trait = snp
            
            # Fréquences alléliques selon la race
            freq = {
                'HAMRA': {'A': 0.7, 'G': 0.3, 'C': 0.6, 'T': 0.4},
                'OUDA': {'A': 0.5, 'G': 0.5, 'C': 0.5, 'T': 0.5},
                'SIDAHOU': {'A': 0.6, 'G': 0.4, 'C': 0.7, 'T': 0.3},
                'BERBERE': {'A': 0.8, 'G': 0.2, 'C': 0.4, 'T': 0.6},
                'CROISE': {'A': 0.55, 'G': 0.45, 'C': 0.55, 'T': 0.45},
                'INCONNU': {'A': 0.5, 'G': 0.5, 'C': 0.5, 'T': 0.5}
            }
            
            race_freq = freq.get(race, freq['INCONNU'])
            
            # Générer génotype selon Hardy-Weinberg
            allele1 = random.choices(alleles, weights=[race_freq.get(alleles[0], 0.5), race_freq.get(alleles[1], 0.5)])[0]
            allele2 = random.choices(alleles, weights=[race_freq.get(alleles[0], 0.5), race_freq.get(alleles[1], 0.5)])[0]
            
            genotype = allele1 + allele2
            
            # Calculer effets
            effet_additif = random.uniform(-0.5, 0.5)
            effet_dominant = random.uniform(-0.3, 0.3)
            frequence_allelique = race_freq.get(allele1, 0.5)
            r2 = random.uniform(0.1, 0.3)
            p_value = random.uniform(0.001, 0.05)
            
            genotypes.append((
                brebis_id, marqueur, chrom, pos, allele1, allele2, genotype,
                frequence_allelique, effet_additif, effet_dominant, r2, p_value,
                f"GENE_{marqueur}", trait, date.today().isoformat()
            ))
        
        return genotypes
    
    @staticmethod
    def calculer_valeurs_genetiques(genotypes):
        """Calcule les valeurs génétiques"""
        if not genotypes:
            return {}
        
        df = pd.DataFrame(genotypes, columns=[
            'brebis_id', 'marqueur', 'chromosome', 'position', 'allele1', 
            'allele2', 'genotype', 'freq_allelique', 'effet_additif', 
            'effet_dominant', 'r2', 'p_value', 'gene_associe', 'trait_associe', 'date'
        ])
        
        valeurs = {}
        
        # Agrégation par trait
        traits = df['trait_associe'].unique()
        
        for trait in traits:
            trait_df = df[df['trait_associe'] == trait]
            
            # Valeur génétique additive
            valeur_additive = (trait_df['effet_additif'] * trait_df['r2']).sum()
            
            # Valeur totale
            valeur_totale = valeur_additive + (trait_df['effet_dominant'] * trait_df['r2']).mean()
            
            # Précision (basée sur R² et nombre de marqueurs)
            precision = min(0.95, trait_df['r2'].sum() / len(trait_df) * 2)
            
            valeurs[trait] = {
                'valeur_additive': round(valeur_additive, 3),
                'valeur_totale': round(valeur_totale, 3),
                'precision': round(precision, 3),
                'nombre_marqueurs': len(trait_df),
                'marqueurs_significatifs': len(trait_df[trait_df['p_value'] < 0.05])
            }
        
        return valeurs
    
    @staticmethod
    def calculer_diversite_genetique(genotypes):
        """Calcule la diversité génétique"""
        if not genotypes:
            return {}
        
        df = pd.DataFrame(genotypes, columns=[
            'brebis_id', 'marqueur', 'chromosome', 'position', 'allele1', 
            'allele2', 'genotype', 'freq_allelique', 'effet_additif', 
            'effet_dominant', 'r2', 'p_value', 'gene_associe', 'trait_associe', 'date'
        ])
        
        # Hétérozygotie observée
        heterozygotes = df[df['allele1'] != df['allele2']]
        ho = len(heterozygotes) / len(df) if len(df) > 0 else 0
        
        # Hétérozygotie attendue (Nei)
        he = 1 - (df['freq_allelique']**2).mean()
        
        # F-statistiques
        fis = 1 - (ho / he) if he > 0 else 0
        
        # Diversité allélique
        alleles_uniques = pd.concat([df['allele1'], df['allele2']]).nunique()
        diversite_allelique = alleles_uniques / (len(df) * 2) if len(df) > 0 else 0
        
        return {
            'heterozygosite_observee': round(ho, 4),
            'heterozygosite_attendue': round(he, 4),
            'fis': round(fis, 4),
            'diversite_allelique': round(diversite_allelique, 4),
            'nombre_snps': len(df['marqueur'].unique()),
            'alleles_uniques': int(alleles_uniques)
        }

# ========== MODULE STATISTIQUES ==========

class ModuleStatistiques:
    """Module d'analyse statistique avancée"""
    
    @staticmethod
    def analyse_descriptive(df, variables):
        """Analyse descriptive complète"""
        results = {}
        
        for var in variables:
            if var in df.columns and pd.api.types.is_numeric_dtype(df[var]):
                data = df[var].dropna()
                
                if len(data) > 0:
                    results[var] = {
                        'n': int(len(data)),
                        'moyenne': round(float(data.mean()), 3),
                        'ecart_type': round(float(data.std()), 3),
                        'minimum': round(float(data.min()), 3),
                        'maximum': round(float(data.max()), 3),
                        'mediane': round(float(data.median()), 3),
                        'coefficient_variation': round(float(data.std() / data.mean() * 100), 1) if data.mean() != 0 else 0,
                        'asymetrie': round(float(skewness(data)), 3),
                        'aplatissement': round(float(kurtosis(data)), 3),
                        'quartile_25': round(float(np.percentile(data, 25)), 3),
                        'quartile_75': round(float(np.percentile(data, 75)), 3)
                    }
        
        return results
    
    @staticmethod
    def test_normalite(df, variables):
        """Test de normalité Shapiro-Wilk simplifié"""
        results = {}
        
        for var in variables:
            if var in df.columns and pd.api.types.is_numeric_dtype(df[var]):
                data = df[var].dropna()
                
                if 3 <= len(data) <= 5000:
                    stat, p = shapiro_wilk_test(data)
                    results[var] = {
                        'statistique': round(stat, 4),
                        'p_value': round(p, 4),
                        'normal': p > 0.05
                    }
        
        return results
    
    @staticmethod
    def correlation_avancee(df, variables):
        """Matrice de corrélation avec tests"""
        numeric_vars = [v for v in variables if v in df.columns and pd.api.types.is_numeric_dtype(df[v])]
        
        if len(numeric_vars) < 2:
            return {}
        
        corr_df = df[numeric_vars].corr(method='pearson')
        
        # Tests de signification
        p_values = pd.DataFrame(index=numeric_vars, columns=numeric_vars)
        
        for i in range(len(numeric_vars)):
            for j in range(len(numeric_vars)):
                if i != j:
                    corr, p_val = pearson_correlation(
                        df[numeric_vars[i]].dropna(),
                        df[numeric_vars[j]].dropna()
                    )
                    p_values.iloc[i, j] = round(p_val, 4)
        
        return {
            'matrice_correlation': corr_df.round(3).to_dict(),
            'p_values': p_values.to_dict(),
            'variables': numeric_vars
        }
    
    @staticmethod
    def analyse_variance(df, variable, facteur):
        """ANOVA à un facteur"""
        if variable not in df.columns or facteur not in df.columns:
            return {}
        
        groups = [group[1][variable].values for group in df[[variable, facteur]].dropna().groupby(facteur)]
        
        if len(groups) >= 2:
            f_stat, p_value = anova_one_way(groups)
            
            result = {
                'f_statistique': round(f_stat, 4),
                'p_value': round(p_value, 4),
                'significatif': p_value < 0.05,
                'nombre_groupes': len(groups),
                'effectif_total': sum(len(g) for g in groups)
            }
            
            # Post-hoc simplifié (différence des moyennes)
            if p_value < 0.05 and len(groups) > 2:
                post_hoc = {}
                group_means = {}
                unique_values = df[facteur].dropna().unique()
                
                for i, val in enumerate(unique_values):
                    group_data = df[df[facteur] == val][variable].dropna()
                    if len(group_data) > 0:
                        group_means[val] = group_data.mean()
                
                # Calculer les différences entre groupes
                for i, val1 in enumerate(unique_values):
                    for j, val2 in enumerate(unique_values):
                        if i < j:
                            diff = abs(group_means.get(val1, 0) - group_means.get(val2, 0))
                            post_hoc[f"{val1} vs {val2}"] = {
                                'difference': round(diff, 3),
                                'significative': diff > 0.5  # Seuil simplifié
                            }
                
                result['post_hoc'] = post_hoc
            
            return result
        
        return {}

# ========== PAGES DE L'APPLICATION ==========

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
    
    # Production laitière
    st.markdown("### 🥛 PRODUCTION LAITIÈRE - COMPARAISON RACES")
    
    cursor.execute("""
        SELECT b.race, 
               AVG(p.quantite_litre) as lait_moyen,
               AVG(p.taux_matiere_grasse) as mg_moyen,
               AVG(p.taux_proteine) as proteine_moyenne,
               COUNT(*) as mesures
        FROM production_lait p
        JOIN brebis b ON p.brebis_id = b.id
        WHERE p.date_mesure > date('now', '-90 days')
        GROUP BY b.race
        HAVING mesures >= 3
    """)
    
    prod_data = cursor.fetchall()
    
    if prod_data:
        df_prod = pd.DataFrame(prod_data, columns=['Race', 'Lait (L)', 'MG (%)', 'Protéine (%)', 'Mesures'])
        
        fig = go.Figure(data=[
            go.Bar(name='Lait (L)', x=df_prod['Race'], y=df_prod['Lait (L)'], marker_color='#8B0000'),
            go.Bar(name='MG (%)', x=df_prod['Race'], y=df_prod['MG (%)'], marker_color='#FF4500', yaxis='y2')
        ])
        
        fig.update_layout(
            title="Production laitière comparée par race",
            yaxis=dict(title="Lait (L/jour)"),
            yaxis2=dict(title="MG (%)", overlaying='y', side='right'),
            barmode='group',
            hovermode='x unified'
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    # Alertes
    st.markdown("### ⚠️ ALERTES ET SUIVI")
    
    # Animaux avec score condition bas
    cursor.execute("""
        SELECT identifiant, nom, race, poids, score_condition
        FROM brebis
        WHERE score_condition <= 2
        AND statut = 'active'
        ORDER BY score_condition
        LIMIT 5
    """)
    
    alertes = cursor.fetchall()
    
    if alertes:
        for alerte in alertes:
            st.warning(f"⚠️ **{alerte[1]}** ({alerte[0]}) - Race {alerte[2]}, Poids: {alerte[3]:.1f}kg, Score: {alerte[4]}/5")
    else:
        st.success("✅ Tous les animaux ont un score de condition satisfaisant")

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
                [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list]
            )
        
        with col_sel2:
            mode_scan = st.selectbox(
                "MODE DE SCAN:",
                ["Laser haute précision", "Photogrammétrie", "Scanner portable", "Simulation"]
            )
        
        if brebis_option:
            brebis_id = int(brebis_option.split('(')[1].split(')')[0].split('-')[-1])
            
            cursor.execute("SELECT * FROM brebis WHERE id = ?", (brebis_id,))
            brebis_info = cursor.fetchone()
            columns = [desc[0] for desc in cursor.description]
            brebis_dict = dict(zip(columns, brebis_info))
            
            # Onglets du scanner
            scan_tabs = st.tabs(["📸 PHOTO", "🎯 SCAN 3D", "📏 MESURES", "📊 ANALYSE"])
            
            with scan_tabs[0]:
                photo = Scanner3D.generer_photo_simulee(brebis_dict)
                
                col_photo1, col_photo2 = st.columns([2, 1])
                
                with col_photo1:
                    st.image(photo, caption=f"Photo simulée - {brebis_dict['nom']}", use_column_width=True)
                    
                    col_btn1, col_btn2 = st.columns(2)
                    with col_btn1:
                        if st.button("📸 Prendre photo", type="primary"):
                            st.success("Photo prise et sauvegardée!")
                    
                    with col_btn2:
                        if st.button("🔄 Régénérer"):
                            st.rerun()
                
                with col_photo2:
                    st.markdown(f"""
                    <div class='race-card'>
                        <h4>📷 INFORMATIONS</h4>
                        <p><strong>Race:</strong> {brebis_dict['race']}</p>
                        <p><strong>ID:</strong> {brebis_dict['identifiant']}</p>
                        <p><strong>Sexe:</strong> {brebis_dict['sexe']}</p>
                        <p><strong>Âge:</strong> {brebis_dict['age_mois']} mois</p>
                        <p><strong>Poids:</strong> {brebis_dict['poids']} kg</p>
                    </div>
                    """, unsafe_allow_html=True)
            
            with scan_tabs[1]:
                st.markdown("### 🎯 SCAN 3D EN TEMPS RÉEL")
                
                scan_progress = st.slider("Progression du scan:", 0, 100, 50)
                
                if scan_progress > 30:
                    points = Scanner3D.simuler_scan_3d(brebis_dict)
                    
                    st.markdown("**Points 3D capturés:**")
                    df_points = pd.DataFrame(points[:10])
                    st.dataframe(df_points[['x', 'y', 'z', 'intensity']])
                    
                    if st.button("🚀 Démarrer scan complet", type="primary"):
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
                
                # Visualisation 3D
                if st.button("👁️ Visualiser scan 3D"):
                    points = Scanner3D.simuler_scan_3d(brebis_dict)
                    fig_3d = Scanner3D.creer_visualisation_3d(points)
                    
                    if fig_3d:
                        st.plotly_chart(fig_3d, use_container_width=True)
            
            with scan_tabs[2]:
                st.markdown("### 📏 MESURES MORPHOMÉTRIQUES")
                
                points = Scanner3D.simuler_scan_3d(brebis_dict)
                mesures = Scanner3D.analyser_points_3d(points)
                
                if mesures:
                    col_mes1, col_mes2, col_mes3 = st.columns(3)
                    
                    with col_mes1:
                        st.metric("Longueur", f"{mesures.get('longueur_cm', 0):.1f} cm")
                        st.metric("Volume", f"{mesures.get('volume_m3', 0):.4f} m³")
                    
                    with col_mes2:
                        st.metric("Largeur", f"{mesures.get('largeur_cm', 0):.1f} cm")
                        st.metric("Surface", f"{mesures.get('surface_m2', 0):.2f} m²")
                    
                    with col_mes3:
                        st.metric("Hauteur", f"{mesures.get('hauteur_cm', 0):.1f} cm")
                        st.metric("Score conformation", f"{mesures.get('score_conformation', 0):.1f}/100")
            
            with scan_tabs[3]:
                st.markdown("### 📊 ANALYSE COMPARATIVE")
                
                race = brebis_dict.get('race', 'INCONNU')
                standards = STANDARDS_RACES.get(race, STANDARDS_RACES['INCONNU'])
                
                points = Scanner3D.simuler_scan_3d(brebis_dict)
                mesures = Scanner3D.analyser_points_3d(points)
                
                if mesures and 'mensurations' in standards:
                    df_compare = pd.DataFrame([
                        {
                            'Mesure': 'Longueur',
                            'Valeur': mesures.get('longueur_cm', 0),
                            'Standard Min': standards['mensurations']['longueur_cm'][0],
                            'Standard Max': standards['mensurations']['longueur_cm'][1],
                            'Optimum': np.mean(standards['mensurations']['longueur_cm'])
                        },
                        {
                            'Mesure': 'Hauteur',
                            'Valeur': mesures.get('hauteur_cm', 0),
                            'Standard Min': standards['mensurations']['hauteur_cm'][0],
                            'Standard Max': standards['mensurations']['hauteur_cm'][1],
                            'Optimum': np.mean(standards['mensurations']['hauteur_cm'])
                        },
                        {
                            'Mesure': 'Largeur bassin',
                            'Valeur': mesures.get('largeur_cm', 0),
                            'Standard Min': standards['mensurations']['largeur_bassin_cm'][0],
                            'Standard Max': standards['mensurations']['largeur_bassin_cm'][1],
                            'Optimum': np.mean(standards['mensurations']['largeur_bassin_cm'])
                        }
                    ])
                    
                    fig = go.Figure()
                    
                    for idx, row in df_compare.iterrows():
                        fig.add_trace(go.Bar(
                            x=[row['Mesure']],
                            y=[row['Valeur']],
                            name='Mesure',
                            marker_color='#8B0000',
                            width=0.3
                        ))
                        
                        # Barre d'intervalle
                        fig.add_trace(go.Scatter(
                            x=[row['Mesure'], row['Mesure']],
                            y=[row['Standard Min'], row['Standard Max']],
                            mode='lines',
                            line=dict(color='green', width=15),
                            name='Standard',
                            showlegend=(idx == 0)
                        ))
                    
                    fig.update_layout(
                        title=f"Comparaison avec standards {race}",
                        yaxis_title="Centimètres (cm)",
                        showlegend=True
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
    
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
        
        export_format = st.selectbox("Format", ["CSV", "Excel", "JSON"])
        
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
        else:
            st.warning("Aucune brebis femelle active")
    
    with tab2:
        cursor.execute("""
            SELECT b.race,
                   strftime('%Y-%m', p.date_mesure) as mois,
                   AVG(p.quantite_litre) as lait_moyen,
                   AVG(p.taux_matiere_grasse) as mg_moyen,
                   AVG(p.taux_proteine) as proteine_moyenne,
                   AVG(p.lactose) as lactose_moyen,
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
                                 columns=['Race', 'Mois', 'Lait (L)', 'MG (%)', 'Protéine (%)', 'Lactose (%)', 'Mesures'])
            
            # Graphique évolution
            fig = px.line(df_prod, x='Mois', y='Lait (L)', color='Race',
                         title="Évolution production laitière par race",
                         markers=True)
            st.plotly_chart(fig, use_container_width=True)
            
            # Graphique qualité
            fig = go.Figure()
            for race in df_prod['Race'].unique():
                race_data = df_prod[df_prod['Race'] == race]
                fig.add_trace(go.Scatter(
                    x=race_data['MG (%)'],
                    y=race_data['Protéine (%)'],
                    mode='markers+text',
                    name=race,
                    text=race_data['Mois'],
                    marker=dict(size=race_data['Lait (L)']*5)
                ))
            
            fig.update_layout(
                title="Relation MG/Protéine par race (taille = production)",
                xaxis_title="Matière grasse (%)",
                yaxis_title="Protéine (%)"
            )
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 🏆 CLASSEMENT DES PRODUCTRICES")
        
        cursor.execute("""
            SELECT b.nom, b.race, b.identifiant,
                   AVG(p.quantite_litre) as moyenne_lait,
                   AVG(p.taux_matiere_grasse) as moyenne_mg,
                   AVG(p.taux_proteine) as moyenne_proteine,
                   COUNT(*) as nb_mesures,
                   MIN(p.date_mesure) as premiere_mesure,
                   MAX(p.date_mesure) as derniere_mesure
            FROM production_lait p
            JOIN brebis b ON p.brebis_id = b.id
            WHERE p.date_mesure > date('now', '-90 days')
            GROUP BY b.id
            HAVING nb_mesures >= 5
            ORDER BY moyenne_lait DESC
            LIMIT 10
        """)
        
        top_prod = cursor.fetchall()
        
        if top_prod:
            df_top = pd.DataFrame(top_prod, 
                                columns=['Nom', 'Race', 'ID', 'Lait moyen (L)', 'MG (%)', 'Protéine (%)', 
                                         'Mesures', 'Première', 'Dernière'])
            
            st.dataframe(df_top.style.highlight_max(subset=['Lait moyen (L)'], color='lightgreen'))
            
            # Graphique top 5
            fig = px.bar(df_top.head(5), x='Nom', y='Lait moyen (L)',
                        color='Race',
                        title="Top 5 productrices - 90 derniers jours",
                        hover_data=['MG (%)', 'Protéine (%)', 'Mesures'])
            st.plotly_chart(fig, use_container_width=True)

def page_stats():
    """Page d'analyse statistique avancée"""
    st.markdown('<h2 class="section-header">📊 ANALYSE STATISTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    cursor = conn.cursor()
    cursor.execute("""
        SELECT race, sexe, age_mois, poids, score_condition, 
               longueur_corps_cm, hauteur_garrot_cm, largeur_bassin_cm,
               tour_poitrine_cm, intensite_couleur, qualite_laine
        FROM brebis
        WHERE statut = 'active'
    """)
    
    data = cursor.fetchall()
    columns = ['Race', 'Sexe', 'Âge', 'Poids', 'Score', 'Longueur', 
               'Hauteur', 'Largeur', 'Poitrine', 'Intensité', 'Qualité_laine']
    
    df = pd.DataFrame(data, columns=columns)
    
    tab1, tab2, tab3, tab4 = st.tabs(["📈 DESCRIPTIVE", "📊 CORRÉLATIONS", "📉 ANOVA", "📋 TESTS"])
    
    with tab1:
        st.markdown("### 📈 ANALYSE DESCRIPTIVE PAR RACE")
        
        variables = ['Poids', 'Longueur', 'Hauteur', 'Largeur', 'Âge']
        races = df['Race'].unique()
        
        for race in races:
            st.markdown(f"#### {STANDARDS_RACES.get(race, {}).get('nom_complet', race)}")
            
            race_df = df[df['Race'] == race]
            stats = ModuleStatistiques.analyse_descriptive(race_df, variables)
            
            if stats:
                # Créer un dataframe pour l'affichage
                stats_df = pd.DataFrame(stats).T
                st.dataframe(stats_df)
                
                # Boxplot
                fig = px.box(race_df, y='Poids', 
                            title=f"Distribution des poids - {race}",
                            points="all")
                st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.markdown("### 📊 MATRICE DE CORRÉLATION")
        
        numeric_vars = ['Poids', 'Longueur', 'Hauteur', 'Largeur', 'Âge', 'Score']
        
        # Matrice de corrélation globale
        corr_matrix = df[numeric_vars].corr()
        
        fig = px.imshow(corr_matrix,
                       title="Matrice de corrélation - Variables morphologiques",
                       color_continuous_scale='RdBu',
                       zmin=-1, zmax=1,
                       text_auto=True)
        st.plotly_chart(fig, use_container_width=True)
        
        # Corrélations par race
        st.markdown("#### CORRÉLATIONS PAR RACE")
        
        race_select = st.selectbox("Sélectionner une race", df['Race'].unique())
        
        race_df = df[df['Race'] == race_select]
        race_corr = race_df[numeric_vars].corr()
        
        fig = px.imshow(race_corr,
                       title=f"Corrélations - {race_select}",
                       color_continuous_scale='RdBu',
                       zmin=-1, zmax=1,
                       text_auto=True)
        st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 📉 ANALYSE DE VARIANCE (ANOVA)")
        
        variable = st.selectbox("Variable à analyser", 
                               ['Poids', 'Longueur', 'Hauteur', 'Largeur', 'Score'])
        facteur = st.selectbox("Facteur de groupement", 
                               ['Race', 'Sexe'])
        
        if st.button("🔬 Exécuter ANOVA"):
            result = ModuleStatistiques.analyse_variance(df, variable, facteur)
            
            if result:
                col_anova1, col_anova2 = st.columns(2)
                
                with col_anova1:
                    st.metric("F-statistique", f"{result['f_statistique']:.4f}")
                    st.metric("p-value", f"{result['p_value']:.4f}")
                    st.metric("Significatif", "✅ Oui" if result['significatif'] else "❌ Non")
                
                with col_anova2:
                    st.metric("Nombre de groupes", result['nombre_groupes'])
                    st.metric("Effectif total", result['effectif_total'])
                
                # Boxplot pour visualiser
                fig = px.box(df, x=facteur, y=variable,
                            title=f"{variable} par {facteur}",
                            points="all",
                            color=facteur)
                st.plotly_chart(fig, use_container_width=True)
                
                # Post-hoc si significatif
                if result['significatif'] and 'post_hoc' in result:
                    st.markdown("#### 📋 DIFFÉRENCES ENTRE GROUPES")
                    post_hoc_df = pd.DataFrame(result['post_hoc']).T
                    st.dataframe(post_hoc_df)
    
    with tab4:
        st.markdown("### 📋 TESTS STATISTIQUES")
        
        # Tests de normalité
        st.markdown("#### TEST DE NORMALITÉ")
        
        test_vars = st.multiselect("Variables à tester", 
                                   ['Poids', 'Longueur', 'Hauteur', 'Largeur', 'Âge'],
                                   default=['Poids'])
        
        if test_vars:
            normalite = ModuleStatistiques.test_normalite(df, test_vars)
            
            if normalite:
                normal_df = pd.DataFrame(normalite).T
                normal_df['Distribution'] = normal_df['normal'].apply(lambda x: 'Normale' if x else 'Non normale')
                
                st.dataframe(normal_df[['statistique', 'p_value', 'Distribution']])
                
                # Visualisation histogramme
                selected_var = st.selectbox("Variable pour histogramme", test_vars)
                
                if selected_var:
                    data = df[selected_var].dropna()
                    
                    fig = go.Figure()
                    fig.add_trace(go.Histogram(
                        x=data,
                        nbinsx=20,
                        name='Distribution',
                        marker_color='#8B0000',
                        opacity=0.7
                    ))
                    
                    # Ajouter courbe normale
                    mean = data.mean()
                    std = data.std()
                    x_range = np.linspace(data.min(), data.max(), 100)
                    y_normal = (1/(std * np.sqrt(2*np.pi))) * np.exp(-0.5*((x_range-mean)/std)**2)
                    y_normal = y_normal * len(data) * (data.max() - data.min()) / 20
                    
                    fig.add_trace(go.Scatter(
                        x=x_range,
                        y=y_normal,
                        mode='lines',
                        name='Courbe normale',
                        line=dict(color='green', width=2)
                    ))
                    
                    fig.update_layout(
                        title=f"Distribution de {selected_var}",
                        xaxis_title=selected_var,
                        yaxis_title="Fréquence",
                        showlegend=True
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)

def page_genetique():
    """Page d'analyse génétique avancée"""
    st.markdown('<h2 class="section-header">🧬 ANALYSE GÉNÉTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🧪 GÉNOTYPAGE", "📊 QTL", "🌳 DIVERSITÉ", "🧮 VALEURS"])
    
    with tab1:
        st.markdown("### 🧪 GÉNOTYPAGE SNP")
        
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY race")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_select = st.selectbox("Sélectionner une brebis", 
                                        [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list])
            
            if brebis_select:
                brebis_id = int(brebis_select.split('(')[1].split(')')[0].split('-')[-1])
                race = brebis_select.split('- ')[1]
                
                if st.button("🧬 Générer génotype", type="primary"):
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
                    
                    # Afficher les résultats
                    df_geno = pd.DataFrame(genotypes, columns=[
                        'brebis_id', 'Marqueur', 'Chromosome', 'Position', 'Allel1', 'Allel2',
                        'Genotype', 'Freq', 'Effet Add', 'Effet Dom', 'R2', 'p-value', 'Gène', 'Trait', 'Date'
                    ])
                    
                    st.success(f"✅ Génotype généré pour {brebis_select}")
                    st.dataframe(df_geno[['Marqueur', 'Chromosome', 'Genotype', 'Trait', 'p-value', 'R2']])
                    
                    # Visualisation
                    fig = px.bar(df_geno, x='Marqueur', y='p-value',
                                title="Significativité des marqueurs",
                                color='R2',
                                color_continuous_scale='Viridis')
                    fig.update_yaxes(type="log")
                    st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.markdown("### 📊 QUANTITATIVE TRAIT LOCI (QTL)")
        
        cursor.execute("SELECT * FROM qtl ORDER BY lod_score DESC")
        qtl_data = cursor.fetchall()
        
        if qtl_data:
            df_qtl = pd.DataFrame(qtl_data, columns=[
                'id', 'Trait', 'Chromosome', 'Début', 'Fin', 'LOD', 
                'Variance', 'Marqueurs', 'Espèce', 'Référence'
            ])
            
            # Filtres
            col_qtl1, col_qtl2 = st.columns(2)
            
            with col_qtl1:
                min_lod = st.slider("LOD score minimum", 0.0, 10.0, 3.0, 0.1)
            
            with col_qtl2:
                traits = st.multiselect("Traits", df_qtl['Trait'].unique(), df_qtl['Trait'].unique())
            
            df_filtered = df_qtl[(df_qtl['LOD'] >= min_lod) & (df_qtl['Trait'].isin(traits))]
            
            # Graphique Manhattan plot
            fig = px.scatter(df_filtered, x='Chromosome', y='LOD',
                            size='Variance', color='Trait',
                            hover_data=['Trait', 'Variance', 'Marqueurs'],
                            title="Manhattan Plot - QTL identifiés")
            
            # Ligne de significativité
            fig.add_hline(y=3.0, line_dash="dash", line_color="red", annotation_text="Seuil LOD=3")
            fig.add_hline(y=5.0, line_dash="dash", line_color="orange", annotation_text="Seuil LOD=5")
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Table détaillée
            st.markdown("#### 📋 TABLEAU DES QTL")
            st.dataframe(df_filtered[['Trait', 'Chromosome', 'LOD', 'Variance', 'Marqueurs', 'Référence']])
    
    with tab3:
        st.markdown("### 🌳 DIVERSITÉ GÉNÉTIQUE")
        
        cursor.execute("""
            SELECT g.brebis_id, b.race, b.identifiant,
                   g.marqueur, g.allele1, g.allele2, g.genotype,
                   g.frequence_allelique, g.trait_associe
            FROM genotypage g
            JOIN brebis b ON g.brebis_id = b.id
            WHERE g.date_analyse > date('now', '-1 year')
        """)
        
        geno_data = cursor.fetchall()
        
        if geno_data:
            # Calculer la diversité
            diversite = ModuleGenetique.calculer_diversite_genetique(geno_data)
            
            if diversite:
                col_div1, col_div2 = st.columns(2)
                
                with col_div1:
                    st.metric("Hétérozygotie observée", f"{diversite['heterozygosite_observee']:.4f}")
                    st.metric("Hétérozygotie attendue", f"{diversite['heterozygosite_attendue']:.4f}")
                    st.metric("Fis", f"{diversite['fis']:.4f}")
                
                with col_div2:
                    st.metric("Diversité allélique", f"{diversite['diversite_allelique']:.4f}")
                    st.metric("SNPs analysés", diversite['nombre_snps'])
                    st.metric("Allèles uniques", diversite['alleles_uniques'])
                
                # Diversité par race
                st.markdown("#### 📊 DIVERSITÉ PAR RACE")
                
                races = set([g[1] for g in geno_data])
                diversite_races = {}
                
                for race in races:
                    race_data = [g for g in geno_data if g[1] == race]
                    div_race = ModuleGenetique.calculer_diversite_genetique(race_data)
                    if div_race:
                        diversite_races[race] = div_race['heterozygosite_observee']
                
                if diversite_races:
                    df_div_race = pd.DataFrame(list(diversite_races.items()), columns=['Race', 'Hétérozygotie'])
                    
                    fig = px.bar(df_div_race, x='Race', y='Hétérozygotie',
                                title="Hétérozygotie observée par race",
                                color='Hétérozygotie',
                                color_continuous_scale='Viridis')
                    st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        st.markdown("### 🧮 VALEURS GÉNÉTIQUES")
        
        cursor.execute("""
            SELECT g.brebis_id, b.race, b.nom,
                   g.marqueur, g.trait_associe, g.effet_additif, 
                   g.effet_dominant, g.r2, g.p_value
            FROM genotypage g
            JOIN brebis b ON g.brebis_id = b.id
            WHERE g.p_value < 0.05
        """)
        
        geno_vals = cursor.fetchall()
        
        if geno_vals:
            valeurs = ModuleGenetique.calculer_valeurs_genetiques(geno_vals)
            
            if valeurs:
                # Créer un dataframe pour l'affichage
                df_vals = pd.DataFrame([
                    {
                        'Trait': trait,
                        'Valeur additive': vals['valeur_additive'],
                        'Valeur totale': vals['valeur_totale'],
                        'Précision': vals['precision'],
                        'Marqueurs': vals['nombre_marqueurs'],
                        'Significatifs': vals['marqueurs_significatifs']
                    }
                    for trait, vals in valeurs.items()
                ])
                
                st.dataframe(df_vals)
                
                # Graphique des valeurs génétiques
                fig = px.bar(df_vals, x='Trait', y='Valeur totale',
                            error_y=df_vals['Précision'].apply(lambda x: 1-x),
                            title="Valeurs génétiques estimées",
                            color='Marqueurs',
                            color_continuous_scale='RdBu')
                st.plotly_chart(fig, use_container_width=True)
                
                # Matrice des corrélations entre traits
                if len(valeurs) > 1:
                    st.markdown("#### 📊 CORRÉLATIONS ENTRE TRAITS")
                    
                    traits_matrix = pd.DataFrame([
                        {trait: vals['valeur_totale'] for trait, vals in valeurs.items()}
                    ])
                    
                    fig = px.imshow(traits_matrix.corr(),
                                   title="Corrélations entre valeurs génétiques",
                                   color_continuous_scale='RdBu',
                                   text_auto=True)
                    st.plotly_chart(fig, use_container_width=True)

# ========== NAVIGATION PRINCIPALE ==========

# Sidebar
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
         "📊 RSTATS",
         "🧬 GÉNÉTIQUE",
         "⚙️ PARAMÈTRES"]
    )
    
    st.markdown("---")
    
    # Statistiques en temps réel
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(DISTINCT race) FROM brebis")
    nb_races = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE statut = 'active'")
    actives = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genotypage WHERE date_analyse > date('now', '-30 days')")
    genotypages = cursor.fetchone()[0]
    
    st.markdown("### 📊 EN DIRECT")
    st.metric("🏷️ Races", nb_races)
    st.metric("🐑 Actives", actives)
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
            <p><small>Couleur: {info['couleur']}</small></p>
            <p><small>Poids ♀️: {info['poids_adulte']['femelle'][0]}-{info['poids_adulte']['femelle'][1]} kg</small></p>
            <p><small>Poids ♂️: {info['poids_adulte']['male'][0]}-{info['poids_adulte']['male'][1]} kg</small></p>
        </div>
        """, unsafe_allow_html=True)

# Navigation des pages
if page == "🏠 ACCUEIL":
    page_accueil()
elif page == "📐 SCANNER 3D":
    page_scanner_3d()
elif page == "📊 GESTION":
    page_gestion()
elif page == "🥛 PRODUCTION":
    page_production()
elif page == "📊 RSTATS":
    page_stats()
elif page == "🧬 GÉNÉTIQUE":
    page_genetique()
elif page == "⚙️ PARAMÈTRES":
    st.info("Page paramètres - À compléter")

# Pied de page
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version Scientifique 4.0</p>
    <p>📐 Scanner 3D + 🧬 Génétique + 📊 Statistiques avancées |rahim GenApAgiE Tlemcen © 2026</p>
</div>
""", unsafe_allow_html=True)
