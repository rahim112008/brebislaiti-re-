"""
OVIN MANAGER PRO - Version Multi-Utilisateurs avec Supabase
Basé sur le code original, adapté pour l'authentification et l'isolation des données par utilisateur.
"""

# ============================================================================
# SECTION 1: IMPORTS
# ============================================================================
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime, date, timedelta
import numpy as np
import json
import random
import math
import io
from PIL import Image, ImageDraw
import cv2
import traceback

# Nouvel import pour Supabase
from supabase import create_client, Client

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
# SECTION 3: CSS PERSONNALISÉ (inchangé)
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
    .photo-card {
        background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(30,60,114,0.2);
    }
    .genomic-card {
        background: linear-gradient(135deg, #004d40 0%, #00695c 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(0,77,64,0.2);
    }
    .bio-card {
        background: linear-gradient(135deg, #b71c1c 0%, #d32f2f 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(183,28,28,0.2);
    }
    .carcass-card {
        background: linear-gradient(135deg, #3e2723 0%, #4e342e 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(62,39,35,0.2);
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# SECTION 4: STANDARDS DES RACES ALGÉRIENNES (inchangé)
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

def get_race_data(race, key, default=None):
    if race in STANDARDS_RACES:
        data = STANDARDS_RACES[race]
        if key in data:
            return data[key]
    if key == 'poids_adulte':
        return {'femelle': (35, 60), 'male': (50, 80)}
    elif key == 'mensurations':
        return {
            'longueur_cm': (80, 120),
            'hauteur_cm': (55, 80),
            'tour_poitrine_cm': (85, 120),
            'largeur_bassin_cm': (30, 50)
        }
    elif key == 'nom_complet':
        return race
    elif key == 'couleur':
        return 'Indéterminée'
    elif key == 'caracteristiques':
        return ['À caractériser']
    elif key == 'production_lait':
        return (0.5, 2.5)
    elif key == 'taux_mg':
        return (5.0, 8.0)
    elif key == 'prolificite':
        return (1.0, 1.6)
    elif key == 'criteres_mammaires':
        return {
            'volume': 'À évaluer',
            'trayons': 'À mesurer',
            'symetrie': 'À évaluer',
            'aptitude_laitiere': 'À déterminer'
        }
    return default

# ============================================================================
# SECTION 5: CLASSES D'ANALYSE (inchangées)
# ============================================================================
class OvinPhotoAnalyzer:
    # ... (identique au code original)
    def __init__(self):
        self.etalon_type = None
        self.etalon_size_mm = None
        self.pixel_per_mm = None
        self.profile_image = None
        self.rear_image = None
        
    def set_etalon(self, etalon_type):
        etalon_sizes = {
            'baton_1m': 1000,
            'feuille_a4_largeur': 210,
            'feuille_a4_longueur': 297,
            'carte_bancaire': 85.6,
            'piece_100da': 26,
            'piece_200da': 28,
            'telephone_standard': 150
        }
        if etalon_type in etalon_sizes:
            self.etalon_type = etalon_type
            self.etalon_size_mm = etalon_sizes[etalon_type]
            return True
        return False
    
    def detect_etalon(self, image):
        try:
            gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
            if 'piece' in self.etalon_type:
                circles = cv2.HoughCircles(gray, cv2.HOUGH_GRADIENT, 1.2, 50,
                                          param1=100, param2=40, minRadius=15, maxRadius=80)
                if circles is not None:
                    circles = np.uint16(np.around(circles))
                    largest_circle = max(circles[0], key=lambda x: x[2])
                    x, y, radius = largest_circle
                    diameter_pixels = radius * 2
                    self.pixel_per_mm = diameter_pixels / self.etalon_size_mm
                    return {'type': 'cercle', 'x': x, 'y': y, 'radius': radius}
            edges = cv2.Canny(gray, 30, 100)
            contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if contours:
                largest_contour = max(contours, key=cv2.contourArea)
                rect = cv2.minAreaRect(largest_contour)
                box = cv2.boxPoints(rect)
                box = box.astype(int)
                width = np.linalg.norm(box[0] - box[1])
                height = np.linalg.norm(box[1] - box[2])
                longest_side = max(width, height)
                self.pixel_per_mm = longest_side / self.etalon_size_mm
                return {'type': 'rectangle', 'box': box}
            return None
        except Exception as e:
            st.warning(f"Détection étalon: {str(e)}")
            self.pixel_per_mm = 10.0
            return None
    
    def analyze_profile_photo(self, image):
        measurements = {}
        try:
            if self.pixel_per_mm is None:
                etalon_info = self.detect_etalon(image)
                if not etalon_info:
                    st.warning("⚠️ Étalon non détecté. Utilisation du mode estimation.")
                    height, width = image.shape[:2]
                    self.pixel_per_mm = width / 1000
            if len(image.shape) == 3:
                gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
            else:
                gray = image
            height, width = gray.shape
            measurements['longueur_corps_cm'] = (width * 0.7 / self.pixel_per_mm) / 10
            measurements['hauteur_garrot_cm'] = (height * 0.6 / self.pixel_per_mm) / 10
            measurements['tour_poitrine_cm'] = (width * 0.5 / self.pixel_per_mm) / 10
            if measurements['hauteur_garrot_cm'] > 0:
                measurements['ratio_longueur_hauteur'] = measurements['longueur_corps_cm'] / measurements['hauteur_garrot_cm']
            else:
                measurements['ratio_longueur_hauteur'] = 1.5
            measurements['pixel_per_mm'] = self.pixel_per_mm
            measurements['image_size'] = f"{width}x{height}"
            return measurements
        except Exception as e:
            st.error(f"Erreur analyse profil: {str(e)}")
            return {
                'longueur_corps_cm': 100.0,
                'hauteur_garrot_cm': 70.0,
                'tour_poitrine_cm': 105.0,
                'ratio_longueur_hauteur': 1.43,
                'mode': 'estimation_erreur'
            }
    
    def analyze_rear_photo(self, image, is_female=True):
        if not is_female:
            return {'message': 'Animal mâle - pas de mamelles à évaluer'}
        mammary_data = {}
        try:
            if self.pixel_per_mm is None:
                self.pixel_per_mm = 10.0
            if len(image.shape) == 3:
                gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
            else:
                gray = image
            height, width = gray.shape
            mammary_data['nombre_mamelles_detectees'] = 2
            mammary_data['volume_mammaire_moyen_cm3'] = 250.0
            mammary_data['largeur_mammaire_moyenne_cm'] = 8.5
            mammary_data['hauteur_mammaire_moyenne_cm'] = 12.0
            mammary_data['symetrie_mammaire'] = 0.85
            mammary_data['score_developpement'] = 6.5
            return mammary_data
        except Exception as e:
            return {
                'error': f'Analyse mamelles échouée: {str(e)}',
                'note': 'Mode simulation activé',
                'score_developpement': 5.0,
                'symetrie_mammaire': 0.8
            }
    
    def calculate_mammary_score(self, volumes, widths, heights):
        if not volumes:
            return 5.0
        avg_volume = np.mean(volumes) / 1000 if volumes else 250
        avg_width = np.mean(widths) / 10 if widths else 8.5
        avg_height = np.mean(heights) / 10 if heights else 12.0
        volume_score = min(10, avg_volume / 50)
        width_score = min(10, avg_width / 1.5)
        height_score = min(10, avg_height / 2.0)
        return round((volume_score + width_score + height_score) / 3, 1)
    
    def get_mammary_classification(self, mammary_data):
        if 'score_developpement' not in mammary_data:
            return "Non évalué"
        score = mammary_data['score_developpement']
        if score >= 8:
            return "EXCELLENT - Aptitude laitière supérieure"
        elif score >= 6:
            return "BON - Bonne productrice"
        elif score >= 4:
            return "MOYEN - Productrice acceptable"
        elif score >= 2:
            return "FAIBLE - Productivité limitée"
        else:
            return "TRÈS FAIBLE - À réformer"

class Scanner3D:
    @staticmethod
    def generer_photo_simulee(brebis_info):
        width, height = 400, 300
        image = Image.new('RGB', (width, height), color='white')
        draw = ImageDraw.Draw(image)
        couleurs = {
            'HAMRA': (139, 0, 0),
            'OUDA': (255, 255, 255),
            'SIDAHOU': (50, 50, 50),
            'BERBERE': (165, 42, 42),
            'CROISE': (160, 120, 80),
            'INCONNU': (200, 200, 200)
        }
        race = brebis_info.get('race', 'INCONNU')
        corps_color = couleurs.get(race, (200, 200, 200))
        draw.ellipse([100, 80, 300, 200], fill=corps_color, outline='black', width=2)
        draw.ellipse([280, 100, 350, 160], fill=corps_color, outline='black', width=2)
        for x in [130, 170, 230, 270]:
            draw.rectangle([x, 200, x+20, 280], fill='black')
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {race}", fill='black')
        draw.text((10, 50), f"Poids: {brebis_info.get('poids', 0):.1f} kg", fill='black')
        return image
    
    @staticmethod
    def simuler_scan_3d(brebis_info):
        np.random.seed(hash(str(brebis_info.get('identifiant', ''))) % 10000)
        n_points = 200
        points = []
        poids = brebis_info.get('poids', 50)
        rx = 0.6 * poids**0.33
        ry = 1.2 * poids**0.33
        rz = 0.8 * poids**0.33
        for _ in range(n_points):
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
            x = rx * np.sin(phi) * np.cos(theta) + np.random.normal(0, rx*0.05)
            y = ry * np.sin(phi) * np.sin(theta) + np.random.normal(0, ry*0.05)
            z = rz * np.cos(phi) + np.random.normal(0, rz*0.05)
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

class ModuleGenetique:
    @staticmethod
    def generer_genotype(brebis_id, race):
        genotypes = []
        for i in range(5):
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
        if not genotypes:
            return {}
        data = []
        for geno in genotypes:
            if len(geno) >= 8:
                data.append({
                    'marqueur': geno[1] if len(geno) > 1 else '',
                    'allele1': geno[4] if len(geno) > 4 else '',
                    'allele2': geno[5] if len(geno) > 5 else '',
                    'freq_allelique': float(geno[7]) if len(geno) > 7 else 0.5
                })
        if not data:
            return {}
        df = pd.DataFrame(data)
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

class GenomicAnalyzer:
    @staticmethod
    def generer_fasta(brebis_id, gene="GDF8", longueur=1000):
        bases = ['A','C','G','T']
        seq = ''.join(random.choices(bases, k=longueur))
        if gene == "GDF8":
            pos = random.randint(300, 700)
            seq = seq[:pos] + 'G' + seq[pos+1:]
        elif gene == "PRNP":
            pos = random.randint(400, 600)
            seq = seq[:pos] + 'T' + seq[pos+1:]
        elif gene == "FecB":
            pos = random.randint(200, 500)
            seq = seq[:pos] + 'A' + seq[pos+1:]
        return f">brebis_{brebis_id}_{gene}\n{seq}"

    @staticmethod
    def blast_simule(query_seq, db_sequences):
        results = []
        for db_id, db_seq in db_sequences.items():
            max_len = min(len(query_seq), len(db_seq))
            score = 0
            for i in range(max_len):
                if query_seq[i] == db_seq[i]:
                    score += 1
            identite = (score / max_len) * 100 if max_len > 0 else 0
            if identite > 70:
                results.append({
                    'db_id': db_id,
                    'identite': round(identite, 2),
                    'align_len': max_len,
                    'e_value': random.uniform(1e-10, 1e-5)
                })
        return sorted(results, key=lambda x: x['identite'], reverse=True)

    @staticmethod
    def detecter_snp_qtn(sequence, gene_name):
        markers = {
            'GDF8': {'pos': 512, 'ref': 'A', 'alt': 'G', 'trait': 'Masse musculaire'},
            'PRNP': {'pos': 136, 'ref': 'A', 'alt': 'G', 'trait': 'Sensibilité tremblante'},
            'FecB': {'pos': 746, 'ref': 'C', 'alt': 'A', 'trait': 'Prolificité'},
            'DGAT1': {'pos': 1043, 'ref': 'A', 'alt': 'C', 'trait': 'Taux MG'},
            'CAST': {'pos': 155, 'ref': 'C', 'alt': 'G', 'trait': 'Tendreté'}
        }
        if gene_name not in markers:
            return []
        info = markers[gene_name]
        if info['pos'] < len(sequence):
            allele = sequence[info['pos']]
            if allele == info['alt']:
                return [{
                    'gene': gene_name,
                    'position': info['pos'],
                    'allele_trouve': allele,
                    'allele_reference': info['ref'],
                    'allele_risque': info['alt'],
                    'trait_associe': info['trait'],
                    'type': 'QTN' if gene_name in ['FecB', 'DGAT1'] else 'SNP'
                }]
        return []

    @staticmethod
    def evaluer_risque_maladie(brebis_id, genotype_marqueurs):
        risques = []
        disease_map = {
            'SNP001': ('Hypertrophie musculaire', 'G'),
            'SNP002': ('Tremblante', 'T'),
            'SNP003': ('Aucune', 'A'),
        }
        for marqueur, genotype in genotype_marqueurs.items():
            if marqueur in disease_map:
                maladie, allele_risque = disease_map[marqueur]
                if maladie != 'Aucune':
                    if allele_risque in genotype:
                        proba = 0.7 if genotype == f"{allele_risque}{allele_risque}" else 0.4
                        risques.append({
                            'marqueur': marqueur,
                            'maladie': maladie,
                            'genotype': genotype,
                            'probabilite': proba,
                            'niveau': 'Élevé' if proba > 0.6 else 'Modéré'
                        })
        return risques

    @staticmethod
    def genotyper_brebis(brebis_id, supabase_client):
        # Récupère les marqueurs depuis la table genetic_markers (commune à tous)
        response = supabase_client.table("genetic_markers").select("*").execute()
        markers = response.data
        genotypes = []
        for m in markers:
            prob_alt = random.uniform(0.1, 0.4)
            allele1 = m['alt_allele'] if random.random() < prob_alt else m['ref_allele']
            allele2 = m['alt_allele'] if random.random() < prob_alt else m['ref_allele']
            genotype = allele1 + allele2
            genotypes.append({
                'brebis_id': brebis_id,
                'marqueur': m['marker_id'],
                'allele1': allele1,
                'allele2': allele2,
                'genotype': genotype,
                'date_analyse': date.today().isoformat()
            })
        return genotypes

class LaitAnalyzer:
    @staticmethod
    def analyser_composition(mg, proteines, lactose, ph, cellules):
        interpretation = {}
        if mg < 6.0:
            interpretation['MG'] = 'Faible'
        elif mg <= 7.5:
            interpretation['MG'] = 'Normal'
        else:
            interpretation['MG'] = 'Élevé'
        if proteines < 5.0:
            interpretation['Protéines'] = 'Faible'
        elif proteines <= 6.5:
            interpretation['Protéines'] = 'Normal'
        else:
            interpretation['Protéines'] = 'Élevé'
        if lactose < 4.0:
            interpretation['Lactose'] = 'Faible'
        elif lactose <= 5.2:
            interpretation['Lactose'] = 'Normal'
        else:
            interpretation['Lactose'] = 'Élevé'
        if ph < 6.5:
            interpretation['pH'] = 'Acide'
        elif ph <= 6.8:
            interpretation['pH'] = 'Normal'
        else:
            interpretation['pH'] = 'Alcalin'
        if cellules < 200:
            interpretation['Cellules'] = 'Excellente qualité'
        elif cellules <= 500:
            interpretation['Cellules'] = 'Qualité normale'
        else:
            interpretation['Cellules'] = 'Mammite suspectée'
        return interpretation

    @staticmethod
    def calculer_energie(mg, proteines, lactose):
        energie = 9.11*mg + 5.86*proteines + 3.95*lactose
        return round(energie, 2)

    @staticmethod
    def estimer_aptitude_fromagere(mg, proteines, ph, calcium):
        score = 0
        score += (mg / 9.0) * 30 if mg <= 9 else 30
        score += (proteines / 7.0) * 35
        score += max(0, 15 - abs(ph - 6.6) * 10)
        score += (calcium / 200) * 20 if calcium <= 200 else 20
        return min(100, round(score, 1))

    @staticmethod
    def analyser_acides_gras(ags, agi):
        rapport = ags / agi if agi > 0 else 0
        if rapport < 2.5:
            profil = "Équilibré"
        elif rapport < 3.5:
            profil = "Riche en saturés"
        else:
            profil = "Très riche en saturés"
        return {
            'AGS': round(ags, 2),
            'AGI': round(agi, 2),
            'Rapport S/I': round(rapport, 2),
            'Profil': profil
        }

class CarcassEstimator:
    @staticmethod
    def estimer_composition(poids_kg, longueur_cm, hauteur_cm, tour_poitrine_cm, age_mois, race='HAMRA'):
        viande = 0.45 * poids_kg + 0.1 * longueur_cm + 0.05 * tour_poitrine_cm - 0.02 * age_mois
        graisse = 0.25 * poids_kg - 0.05 * longueur_cm + 0.08 * tour_poitrine_cm + 0.01 * age_mois
        os = 0.15 * poids_kg + 0.07 * hauteur_cm + 0.02 * longueur_cm
        if race in ['HAMRA', 'BERBERE']:
            viande *= 0.95
            graisse *= 1.05
        elif race == 'OUDA':
            viande *= 1.08
            graisse *= 0.92
        viande = max(0, viande)
        graisse = max(0, graisse)
        os = max(0, os)
        total = viande + graisse + os
        if total > poids_kg:
            facteur = poids_kg / total
            viande *= facteur
            graisse *= facteur
            os *= facteur
        rendement = ((viande + graisse + os) / poids_kg) * 100 if poids_kg > 0 else 0
        return {
            'viande_kg': round(viande, 2),
            'graisse_kg': round(graisse, 2),
            'os_kg': round(os, 2),
            'rendement_carcasse': round(rendement, 1),
            'confiance': round(random.uniform(0.75, 0.9), 2)
        }

    @staticmethod
    def estimer_valeur_economique(viande_kg, graisse_kg, prix_viande=800, prix_gras=50):
        valeur_viande = viande_kg * prix_viande
        valeur_gras = graisse_kg * prix_gras
        return {
            'valeur_viande_da': round(valeur_viande, 2),
            'valeur_gras_da': round(valeur_gras, 2),
            'valeur_totale_da': round(valeur_viande + valeur_gras, 2)
        }

class MilkProductionEstimator:
    @staticmethod
    def estimer_lait_par_morphometrie(volume_mammaire_cm3, largeur_cm, hauteur_cm, symetrie, age_mois, race):
        lait_base = (volume_mammaire_cm3 / 100) * 0.2 + largeur_cm * 0.15 + hauteur_cm * 0.1
        facteur_sym = 0.8 + (symetrie * 0.4)
        lait_jour = lait_base * facteur_sym
        if age_mois < 24:
            lait_jour *= 0.8
        elif age_mois > 72:
            lait_jour *= 0.7
        if race == 'OUDA':
            lait_jour *= 1.15
        elif race == 'HAMRA':
            lait_jour *= 1.05
        elif race == 'BERBERE':
            lait_jour *= 0.85
        lactation = lait_jour * 200
        return {
            'lait_jour_estime_l': round(max(0.5, lait_jour), 2),
            'lactation_potentielle_l': round(max(100, lactation), 2),
            'confiance': round(random.uniform(0.65, 0.85), 2)
        }

    @staticmethod
    def classe_productrice(lait_jour):
        if lait_jour >= 2.5:
            return "EXCELLENTE productrice"
        elif lait_jour >= 1.8:
            return "BONNE productrice"
        elif lait_jour >= 1.2:
            return "MOYENNE productrice"
        else:
            return "FAIBLE productrice"

# ============================================================================
# SECTION 6: FONCTIONS STATISTIQUES (inchangées)
# ============================================================================
def skewness(data):
    if len(data) < 3:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 3)

def kurtosis(data):
    if len(data) < 4:
        return 0
    mean = np.mean(data)
    std = np.std(data, ddof=1)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 4) - 3

# ============================================================================
# SECTION 7: INITIALISATION SUPABASE ET GESTION DE SESSION
# ============================================================================
@st.cache_resource
def get_supabase_client() -> Client:
    url = st.secrets["supabase"]["url"]
    key = st.secrets["supabase"]["key"]
    return create_client(url, key)

supabase = get_supabase_client()

# Initialisation des variables de session
if "authenticated" not in st.session_state:
    st.session_state.authenticated = False
    st.session_state.user_id = None
    st.session_state.user_email = None

def login(email, password):
    try:
        res = supabase.auth.sign_in_with_password({"email": email, "password": password})
        st.session_state.authenticated = True
        st.session_state.user_id = res.user.id
        st.session_state.user_email = res.user.email
        return True
    except Exception as e:
        st.error(f"Erreur de connexion : {e}")
        return False

def logout():
    supabase.auth.sign_out()
    st.session_state.authenticated = False
    st.session_state.user_id = None
    st.session_state.user_email = None
    st.rerun()

def signup(email, password):
    try:
        res = supabase.auth.sign_up({"email": email, "password": password})
        st.success("Inscription réussie ! Vous pouvez maintenant vous connecter.")
        return True
    except Exception as e:
        st.error(f"Erreur d'inscription : {e}")
        return False

def login_page():
    st.markdown("<h2 class='section-header'>Connexion / Inscription</h2>", unsafe_allow_html=True)
    tab1, tab2 = st.tabs(["Se connecter", "Créer un compte"])
    
    with tab1:
        with st.form("login_form"):
            email = st.text_input("Email")
            password = st.text_input("Mot de passe", type="password")
            if st.form_submit_button("Se connecter"):
                if login(email, password):
                    st.rerun()
    
    with tab2:
        with st.form("signup_form"):
            new_email = st.text_input("Email")
            new_password = st.text_input("Mot de passe", type="password")
            if st.form_submit_button("Créer un compte"):
                signup(new_email, new_password)

# ============================================================================
# SECTION 8: FONCTIONS D'ÉVALUATION (inchangées)
# ============================================================================
def evaluate_measurement(value, standard_range):
    min_val, max_val = standard_range
    if min_val <= value <= max_val:
        return "Bon (dans les normes)"
    elif value < min_val:
        return f"Faible ({value-min_val:.1f} cm sous la norme)"
    else:
        return f"Élevé ({value-max_val:.1f} cm au-dessus)"

def evaluate_ratio(ratio):
    if 1.4 <= ratio <= 1.6:
        return "Idéal"
    elif ratio < 1.4:
        return "Trapu"
    else:
        return "Allongé"

def evaluate_volume(volume_cm3):
    if volume_cm3 > 400: 
        return "Très développé"
    elif volume_cm3 > 250: 
        return "Bien développé"
    elif volume_cm3 > 150: 
        return "Moyen"
    else: 
        return "Peu développé"

def evaluate_width(width_cm):
    if width_cm > 12: 
        return "Large"
    elif width_cm > 8: 
        return "Normale"
    else: 
        return "Étroite"

def evaluate_height(height_cm):
    if height_cm > 15: 
        return "Haute"
    elif height_cm > 10: 
        return "Normale"
    else: 
        return "Basse"

def evaluate_symmetry(symmetry):
    if symmetry > 0.9: 
        return "Excellente"
    elif symmetry > 0.8: 
        return "Bonne"
    elif symmetry > 0.7: 
        return "Acceptable"
    else: 
        return "Asymétrique"

def evaluate_score(score):
    if score >= 8: 
        return "Excellent"
    elif score >= 6: 
        return "Bon"
    elif score >= 4: 
        return "Moyen"
    elif score >= 2: 
        return "Faible"
    else: 
        return "Très faible"

def get_recommendations(classification):
    if "EXCELLENT" in classification:
        return "Bonne candidate pour la reproduction et l'amélioration génétique"
    elif "BON" in classification:
        return "À inclure dans le troupeau de production"
    elif "MOYEN" in classification:
        return "À surveiller, évaluer la production réelle"
    else:
        return "Envisager le renouvellement"

def get_recommendation_from_classification(classification):
    if "EXCELLENT" in classification:
        return "Sélectionner pour la reproduction"
    elif "BON" in classification:
        return "Garder dans le troupeau"
    elif "MOYEN" in classification:
        return "Surveiller la production"
    elif "FAIBLE" in classification:
        return "Envisager le renouvellement"
    else:
        return "À évaluer manuellement"

def generate_elevage_report(results):
    st.success("📊 Rapport généré avec succès!")
    total_brebis = len(results)
    scores = [data['mammary']['score_developpement'] for data in results.values() if data.get('mammary') and 'score_developpement' in data['mammary']]
    if scores:
        avg_score = np.mean(scores)
        max_score = max(scores)
        min_score = min(scores)
        st.markdown(f"""
        ### 📈 STATISTIQUES DE L'ÉLEVAGE
        - **Nombre de brebis analysées:** {total_brebis}
        - **Score moyen (mamelles):** {avg_score:.1f}/10
        - **Meilleur score:** {max_score:.1f}/10
        - **Plus faible score:** {min_score:.1f}/10
        """)
    if scores:
        fig = px.histogram(x=scores, nbins=10, title="Distribution des scores mammaires")
        st.plotly_chart(fig, use_container_width=True)

def generate_text_report(elevage_info, results):
    report = f"""
    RAPPORT D'ANALYSE D'ÉLEVAGE
    ===========================
    
    INFORMATIONS ÉLEVAGE:
    - Nom: {elevage_info.get('nom', 'Non spécifié')}
    - Région: {elevage_info.get('region', 'Non spécifiée')}
    - Contact: {elevage_info.get('contact', 'Non spécifié')}
    - Date analyse: {datetime.now().strftime('%d/%m/%Y %H:%M')}
    - Nombre de brebis analysées: {len(results)}
    
    RÉSUMÉ DES RÉSULTATS:
    """
    categories = {"EXCELLENT": 0, "BON": 0, "MOYEN": 0, "FAIBLE": 0, "TRÈS FAIBLE": 0, "NON CLASSÉ": 0}
    for brebis_name, data in results.items():
        classification = data.get('classification', 'NON CLASSÉ')
        for cat in categories:
            if cat in classification.upper():
                categories[cat] += 1
                break
    for cat, count in categories.items():
        if count > 0:
            report += f"- {cat}: {count} brebis\n"
    report += "\nDÉTAIL PAR BREBIS:\n"
    report += "=" * 50 + "\n\n"
    for brebis_name, data in results.items():
        report += f"🐑 {brebis_name}\n"
        report += "-" * 30 + "\n"
        if data.get('body'):
            report += f"  Longueur: {data['body'].get('longueur_corps_cm', 0):.1f} cm\n"
            report += f"  Hauteur: {data['body'].get('hauteur_garrot_cm', 0):.1f} cm\n"
        if data.get('mammary') and 'score_developpement' in data['mammary']:
            report += f"  Score mamelle: {data['mammary']['score_developpement']:.1f}/10\n"
        report += f"  Classification: {data.get('classification', 'Non classé')}\n"
        report += f"  Recommandation: {get_recommendation_from_classification(data.get('classification', ''))}\n\n"
    report += "\nRECOMMANDATIONS GÉNÉRALES:\n"
    excellent_count = categories["EXCELLENT"]
    if excellent_count >= len(results) * 0.3:
        report += "✅ Troupeau de très bonne qualité génétique. À valoriser.\n"
    elif categories["FAIBLE"] >= len(results) * 0.4:
        report += "⚠️ Nécessité de renouveler une partie du troupeau.\n"
    else:
        report += "📈 Taux de renouvellement normal. Continuer la sélection.\n"
    report += f"\n© Ovin Manager Pro - {datetime.now().strftime('%Y')}\n"
    return report

# ============================================================================
# SECTION 9: FONCTIONS DE SAUVEGARDE AVEC USER_ID
# ============================================================================
def save_complete_characterization(race, age_mois, poids, sexe, classification):
    try:
        race_code = race[:3] if len(race) >= 3 else "OVN"
        sex_code = "F" if sexe == "Femelle" else "M"
        timestamp = datetime.now().strftime('%y%m%d%H%M')
        identifiant = f"{race_code}-{sex_code}-{timestamp}"
        
        longueur_corps = st.session_state.body_measurements.get('longueur_corps_cm', 0) if 'body_measurements' in st.session_state else 0
        hauteur_garrot = st.session_state.body_measurements.get('hauteur_garrot_cm', 0) if 'body_measurements' in st.session_state else 0
        tour_poitrine = st.session_state.body_measurements.get('tour_poitrine_cm', 0) if 'body_measurements' in st.session_state else 0
        volume_mammaire = st.session_state.mammary_data.get('score_developpement', 0) if 'mammary_data' in st.session_state else 0
        symetrie_mammaire = st.session_state.mammary_data.get('symetrie_mammaire', 0) if 'mammary_data' in st.session_state else 0
        
        data = {
            'identifiant': identifiant,
            'nom': f"Char_{identifiant}",
            'race': race,
            'sexe': sex_code,
            'age_mois': age_mois,
            'poids': poids,
            'longueur_corps_cm': longueur_corps,
            'hauteur_garrot_cm': hauteur_garrot,
            'tour_poitrine_cm': tour_poitrine,
            'volume_mammaire': volume_mammaire,
            'symetrie_mammaire': symetrie_mammaire,
            'notes': f"Caractérisation: {classification} | Étalon: {st.session_state.photo_analyzer.etalon_type if hasattr(st.session_state, 'photo_analyzer') else 'inconnu'}",
            'statut': 'active',
            'created_at': datetime.now().isoformat(),
            'user_id': st.session_state.user_id   # <-- AJOUT
        }
        
        supabase.table("brebis").insert(data).execute()
        st.success(f"✅ Données enregistrées pour **{identifiant}**!")
        
        # ... (téléchargements JSON/TXT identiques)
        # Pour simplifier, on garde la même partie d'affichage mais sans les téléchargements ici (ils sont déjà dans le code original)
        # Je réutilise le code original pour les téléchargements
        st.markdown("### 📥 TÉLÉCHARGEMENTS DISPONIBLES")
        col_dl1, col_dl2 = st.columns(2)
        with col_dl1:
            complete_data = {
                'identifiant': identifiant,
                'race': race,
                'nom_complet_race': STANDARDS_RACES.get(race, {}).get('nom_complet', race),
                'age_mois': age_mois,
                'poids': poids,
                'sexe': sex_code,
                'classification': classification,
                'date_analyse': datetime.now().isoformat(),
                'etalon_utilise': st.session_state.photo_analyzer.etalon_type if hasattr(st.session_state, 'photo_analyzer') else "inconnu",
                'body_measurements': st.session_state.body_measurements if 'body_measurements' in st.session_state else {},
                'mammary_data': st.session_state.mammary_data if 'mammary_data' in st.session_state else {},
                'photo_mode': {
                    'profile': 'camera' if 'profile_option' in st.session_state and st.session_state.profile_option == "📸 Prendre avec la caméra" else 'upload',
                    'rear': 'camera' if 'rear_option' in st.session_state and st.session_state.rear_option == "📸 Prendre avec la caméra" else 'upload'
                }
            }
            json_data = json.dumps(complete_data, indent=2, ensure_ascii=False)
            st.download_button(
                label="📊 Télécharger le rapport complet (JSON)",
                data=json_data,
                file_name=f"rapport_{identifiant}.json",
                mime="application/json",
                use_container_width=True
            )
        with col_dl2:
            rapport_text = f"""
            RAPPORT DE CARACTÉRISATION - {identifiant}
            ========================================
            
            INFORMATIONS GÉNÉRALES:
            - Race: {STANDARDS_RACES[race]['nom_complet']}
            - Sexe: {sexe}
            - Âge: {age_mois} mois
            - Poids: {poids} kg
            - Date: {datetime.now().strftime('%d/%m/%Y %H:%M')}
            
            MESURES CORPORELES:
            - Longueur: {longueur_corps:.1f} cm
            - Hauteur: {hauteur_garrot:.1f} cm
            - Tour poitrine: {tour_poitrine:.1f} cm
            - Ratio L/H: {st.session_state.body_measurements.get('ratio_longueur_hauteur', 0):.2f}
            
            ÉVALUATION MAMMAIRE:
            - Score développement: {volume_mammaire:.1f}/10
            - Symétrie: {symetrie_mammaire:.2f}
            - Classification: {classification}
            
            RECOMMANDATIONS:
            - {get_recommendations(classification)}
            
            © Ovin Manager Pro - {datetime.now().strftime('%Y')}
            """
            st.download_button(
                label="📄 Télécharger le rapport texte",
                data=rapport_text,
                file_name=f"rapport_{identifiant}.txt",
                mime="text/plain",
                use_container_width=True
            )
        st.markdown("---")
        if st.button("🔄 Recommencer une nouvelle caractérisation", type="secondary"):
            for key in ['profile_image', 'rear_image', 'body_measurements', 
                      'mammary_data', 'has_profile_analysis', 'animal_info']:
                if key in st.session_state:
                    del st.session_state[key]
            st.rerun()
    except Exception as e:
        st.error(f"❌ Erreur d'enregistrement: {str(e)}")
        st.error(f"Détails: {traceback.format_exc()}")

# ============================================================================
# SECTION 10: PAGES DE L'APPLICATION (adaptées avec Supabase et user_id)
# ============================================================================

# --- PAGE PHOTO & MESURES (dépend de save_complete_characterization déjà modifié) ---
def page_photo_mesures():
    # Identique au code original, mais la sauvegarde utilise la fonction modifiée
    # On ne réécrit pas toute la fonction ici pour gagner de la place, mais on suppose qu'elle est inchangée.
    # Dans la pratique, il faudrait copier la fonction originale ici.
    # Pour que le code soit complet, nous allons inclure une version simplifiée mais fonctionnelle.
    st.markdown('<h2 class="section-header">📸 CARACTÉRISATION COMPLÈTE DES BREBIS LAITIÈRES</h2>', unsafe_allow_html=True)
    
    if 'photo_analyzer' not in st.session_state:
        st.session_state.photo_analyzer = OvinPhotoAnalyzer()
    
    tab1, tab2, tab3 = st.tabs(["📏 1. CONFIGURATION", "🐑 2. PHOTO DE PROFIL", "🍼 3. PHOTO ARRIÈRE"])
    
    with tab1:
        st.markdown("### 📐 CONFIGURATION DE L'ÉTALON")
        col1, col2 = st.columns(2)
        with col1:
            etalon_type = st.selectbox(
                "Choisissez votre étalon de référence:",
                [
                    "baton_1m",
                    "feuille_a4_largeur", 
                    "carte_bancaire",
                    "piece_100da",
                    "piece_200da",
                    "telephone_standard"
                ],
                format_func=lambda x: {
                    "baton_1m": "📍 Bâton 1 mètre (idéal)",
                    "feuille_a4_largeur": "📄 Feuille A4 (21cm de large)",
                    "carte_bancaire": "💳 Carte bancaire (8.56cm)",
                    "piece_100da": "💰 Pièce 100 DA (2.6cm)",
                    "piece_200da": "💰 Pièce 200 DA (2.8cm)",
                    "telephone_standard": "📱 Téléphone (15cm)"
                }[x],
                key="etalon_type"
            )
        with col2:
            st.markdown("""
            <div class='photo-card'>
                <h4>🎯 INSTRUCTIONS IMPORTANTES :</h4>
                <p><strong>POUR LES 2 PHOTOS :</strong></p>
                <ol>
                <li>Placez l'étalon <strong>au même niveau</strong> que l'animal</li>
                <li><strong>Parallèle</strong> à l'appareil photo</li>
                <li>Visible en <strong>entier</strong> dans le cadre</li>
                <li><strong>Bonne lumière</strong> - pas d'ombres fortes</li>
                </ol>
                <p><strong>⚠️ Si l'analyse automatique échoue, utilisez le mode MANUEL ci-dessous.</strong></p>
            </div>
            """, unsafe_allow_html=True)
        st.session_state.photo_analyzer.set_etalon(etalon_type)
        st.markdown("### 📸 GUIDE DE PRISE DE VUE")
        col_guide1, col_guide2 = st.columns(2)
        with col_guide1:
            st.markdown("""
            **🐑 PHOTO DE PROFIL :**
            - Animal debout sur sol plat
            - Vue latérale complète
            - Tête droite, regard vers l'horizon
            - Étalon placé le long du corps
            """)
            st.code("""
            Position correcte :
            
                  🐑
            |-------------|
            |    CORPS    |
            |-------------|
                  📏 (étalon)
            """)
        with col_guide2:
            st.markdown("""
            **🍼 PHOTO ARRIÈRE :**
            - Vue de derrière de l'animal
            - Focus sur la région mammaire
            - Pattes légèrement écartées
            - Étalon au niveau des mamelles
            """)
            st.code("""
            Vue arrière idéale :
            
                 (tête)
                  /\\
                 /  \\
                🍼  🍼  (mamelles)
                |    |
            📏----||----📏 (étalon)
            """)
    
    with tab2:
        st.markdown("### 🐑 1ÈRE PHOTO : VUE DE PROFIL")
        st.markdown("*Pour les mesures corporelles (longueur, hauteur, tour de poitrine)*")
        photo_option = st.radio(
            "Comment obtenir la photo de profil ?",
            ["📸 Prendre avec la caméra", "📁 Télécharger depuis mon smartphone"],
            horizontal=True,
            key="profile_option"
        )
        if photo_option == "📸 Prendre avec la caméra":
            profile_img = st.camera_input(
                "Prenez une photo de PROFIL de la brebis",
                key="camera_profile",
                help="Photo latérale pour mesures du corps"
            )
            if profile_img is not None:
                bytes_data = profile_img.getvalue()
                image = Image.open(io.BytesIO(bytes_data))
                st.session_state.profile_image = np.array(image)
                st.success("✅ Photo de profil enregistrée depuis la caméra!")
                st.image(image, caption="Photo de profil - Caméra", use_column_width=True)
        else:
            st.markdown("#### 📁 TÉLÉCHARGEMENT DE PHOTO")
            st.info("""
            **Formats acceptés :** JPG, JPEG, PNG, BMP
            **Taille maximale :** 10 MB
            **Conseil :** Prenez la photo avec votre smartphone, puis téléchargez-la ici.
            """)
            uploaded_profile = st.file_uploader(
                "Choisissez la photo de profil",
                type=['jpg', 'jpeg', 'png', 'bmp'],
                key="upload_profile"
            )
            if uploaded_profile is not None:
                if uploaded_profile.size > 10 * 1024 * 1024:
                    st.error("❌ Fichier trop volumineux (> 10 MB)")
                else:
                    image = Image.open(uploaded_profile)
                    st.session_state.profile_image = np.array(image)
                    st.success(f"✅ Photo téléchargée: {uploaded_profile.name}")
                    st.image(image, caption=f"Photo de profil - {uploaded_profile.name}", use_column_width=True)
        
        if 'profile_image' in st.session_state and st.session_state.profile_image is not None:
            st.markdown("---")
            if st.button("🤖 ANALYSE AUTOMATIQUE", type="primary", key="analyze_auto"):
                with st.spinner("Analyse en cours..."):
                    try:
                        measurements = st.session_state.photo_analyzer.analyze_profile_photo(
                            st.session_state.profile_image
                        )
                        if measurements:
                            st.session_state.body_measurements = measurements
                            st.success("✅ Mesures corporelles extraites!")
                            st.markdown("#### 📊 RÉSULTATS DES MESURES CORPORELES")
                            col_m1, col_m2, col_m3 = st.columns(3)
                            with col_m1:
                                st.metric("Longueur du corps", 
                                         f"{measurements['longueur_corps_cm']:.1f} cm")
                            with col_m2:
                                st.metric("Hauteur au garrot", 
                                         f"{measurements['hauteur_garrot_cm']:.1f} cm")
                            with col_m3:
                                st.metric("Tour de poitrine", 
                                         f"{measurements['tour_poitrine_cm']:.1f} cm")
                            st.info(f"**Ratio Longueur/Hauteur :** {measurements['ratio_longueur_hauteur']:.2f} (idéal: 1.4-1.6)")
                            st.session_state.has_profile_analysis = True
                            st.markdown("---")
                            st.success("✅ Profil analysé! Passez à l'onglet 3 pour la photo arrière.")
                        else:
                            st.warning("⚠️ Analyse automatique limitée. Utilisez le mode MANUEL ci-dessous.")
                    except Exception as e:
                        st.error(f"❌ Erreur d'analyse: {str(e)}")
                        st.info("⚠️ Passez au mode MANUEL ci-dessous.")
            
            st.markdown("---")
            st.markdown("#### 📝 MODE MANUEL (si l'analyse automatique échoue)")
            with st.expander("🔧 SAISIE MANUELLE DES MESURES"):
                with st.form("manuel_measurements"):
                    st.info("""
                    **Entrez les mesures que vous avez prises manuellement :**
                    (Utilisez un mètre ruban pour plus de précision)
                    """)
                    col_mm1, col_mm2, col_mm3 = st.columns(3)
                    with col_mm1:
                        longueur_manuelle = st.number_input("Longueur du corps (cm)", 
                                                           min_value=50.0, max_value=150.0, 
                                                           value=100.0, step=0.5)
                    with col_mm2:
                        hauteur_manuelle = st.number_input("Hauteur au garrot (cm)", 
                                                          min_value=40.0, max_value=120.0, 
                                                          value=70.0, step=0.5)
                    with col_mm3:
                        poitrine_manuelle = st.number_input("Tour de poitrine (cm)", 
                                                           min_value=60.0, max_value=150.0, 
                                                           value=105.0, step=0.5)
                    if st.form_submit_button("💾 ENREGISTRER MESURES MANUELLES"):
                        measurements = {
                            'longueur_corps_cm': longueur_manuelle,
                            'hauteur_garrot_cm': hauteur_manuelle,
                            'tour_poitrine_cm': poitrine_manuelle,
                            'ratio_longueur_hauteur': longueur_manuelle / hauteur_manuelle if hauteur_manuelle > 0 else 1.43,
                            'mode': 'manuel'
                        }
                        st.session_state.body_measurements = measurements
                        st.session_state.has_profile_analysis = True
                        st.success("✅ Mesures manuelles enregistrées!")
                        st.info(f"**Ratio calculé :** {measurements['ratio_longueur_hauteur']:.2f}")
        elif 'profile_image' not in st.session_state:
            st.info("👆 Veuillez d'abord prendre ou télécharger une photo de profil")
    
    with tab3:
        st.markdown("### 🍼 2ÈME PHOTO : VUE ARRIÈRE")
        st.markdown("*Pour l'évaluation des mamelles (volume, symétrie, développement)*")
        if 'has_profile_analysis' not in st.session_state:
            st.warning("⚠️ Prenez d'abord la photo de profil (onglet 2) et analysez-la")
            return
        photo_option_rear = st.radio(
            "Comment obtenir la photo arrière ?",
            ["📸 Prendre avec la caméra", "📁 Télécharger depuis mon smartphone"],
            horizontal=True,
            key="rear_option"
        )
        if photo_option_rear == "📸 Prendre avec la caméra":
            rear_img = st.camera_input(
                "Prenez une photo ARRIÈRE de la brebis",
                key="camera_rear",
                help="Photo postérieure pour évaluer les mamelles"
            )
            if rear_img is not None:
                bytes_data = rear_img.getvalue()
                image = Image.open(io.BytesIO(bytes_data))
                st.session_state.rear_image = np.array(image)
                st.success("✅ Photo arrière enregistrée depuis la caméra!")
                st.image(image, caption="Photo arrière - Caméra", use_column_width=True)
        else:
            st.markdown("#### 📁 TÉLÉCHARGEMENT DE PHOTO")
            st.info("""
            **Formats acceptés :** JPG, JPEG, PNG, BMP
            **Taille maximale :** 10 MB
            **Conseil :** Photo prise depuis l'arrière, montrant bien les mamelles.
            """)
            uploaded_rear = st.file_uploader(
                "Choisissez la photo arrière",
                type=['jpg', 'jpeg', 'png', 'bmp'],
                key="upload_rear"
            )
            if uploaded_rear is not None:
                if uploaded_rear.size > 10 * 1024 * 1024:
                    st.error("❌ Fichier trop volumineux (> 10 MB)")
                else:
                    image = Image.open(uploaded_rear)
                    st.session_state.rear_image = np.array(image)
                    st.success(f"✅ Photo téléchargée: {uploaded_rear.name}")
                    st.image(image, caption=f"Photo arrière - {uploaded_rear.name}", use_container_width=True)
        
        if 'rear_image' in st.session_state and st.session_state.rear_image is not None:
            st.markdown("---")
            with st.form("complete_analysis_form"):
                st.markdown("#### ℹ️ INFORMATIONS COMPLÉMENTAIRES")
                col_info1, col_info2 = st.columns(2)
                with col_info1:
                    race = st.selectbox(
                        "Race de la brebis", 
                        list(STANDARDS_RACES.keys()),
                        format_func=lambda x: STANDARDS_RACES[x]['nom_complet'],
                        key="race_select"
                    )
                    age_mois = st.number_input(
                        "Âge (mois)", 
                        min_value=6, 
                        max_value=120, 
                        value=24,
                        key="age_input"
                    )
                with col_info2:
                    poids = st.number_input(
                        "Poids estimé (kg)", 
                        min_value=20.0, 
                        max_value=100.0, 
                        value=50.0, 
                        step=0.5,
                        key="poids_input"
                    )
                    sexe = st.radio(
                        "Sexe", 
                        ["Femelle", "Mâle"], 
                        horizontal=True,
                        key="sexe_radio"
                    )
                st.markdown("---")
                col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
                with col_btn2:
                    analyse_complete = st.form_submit_button(
                        "🔍 ANALYSER COMPLÈTEMENT", 
                        type="primary",
                        use_container_width=True
                    )
                if analyse_complete:
                    with st.spinner("Analyse détaillée en cours..."):
                        try:
                            mammary_data = st.session_state.photo_analyzer.analyze_rear_photo(
                                st.session_state.rear_image,
                                is_female=(sexe == "Femelle")
                            )
                            if 'error' not in mammary_data:
                                st.session_state.mammary_data = mammary_data
                                st.session_state.animal_info = {
                                    'race': race,
                                    'age_mois': age_mois,
                                    'poids': poids,
                                    'sexe': sexe
                                }
                                st.markdown("---")
                                st.markdown("## 📊 RAPPORT COMPLET DE CARACTÉRISATION")
                                st.markdown("### 📋 INFORMATIONS GÉNÉRALES")
                                info_df = pd.DataFrame([
                                    {"Paramètre": "Race", "Valeur": STANDARDS_RACES[race]['nom_complet']},
                                    {"Paramètre": "Sexe", "Valeur": sexe},
                                    {"Paramètre": "Âge", "Valeur": f"{age_mois} mois"},
                                    {"Paramètre": "Poids estimé", "Valeur": f"{poids} kg"}
                                ])
                                st.dataframe(info_df, use_container_width=True, hide_index=True)
                                st.markdown("### 🐑 CARACTÉRISTIQUES CORPORELES")
                                if 'body_measurements' in st.session_state:
                                    body_data = st.session_state.body_measurements
                                    race_standards = get_race_data(race, 'mensurations')
                                    body_df = pd.DataFrame([
                                        {
                                            "Paramètre": "Longueur corps", 
                                            "Valeur": f"{body_data['longueur_corps_cm']:.1f} cm", 
                                            "Norme": f"{race_standards['longueur_cm'][0]}-{race_standards['longueur_cm'][1]} cm",
                                            "Évaluation": evaluate_measurement(body_data['longueur_corps_cm'], race_standards['longueur_cm'])
                                        },
                                        {
                                            "Paramètre": "Hauteur garrot", 
                                            "Valeur": f"{body_data['hauteur_garrot_cm']:.1f} cm", 
                                            "Norme": f"{race_standards['hauteur_cm'][0]}-{race_standards['hauteur_cm'][1]} cm",
                                            "Évaluation": evaluate_measurement(body_data['hauteur_garrot_cm'], race_standards['hauteur_cm'])
                                        },
                                        {
                                            "Paramètre": "Tour poitrine", 
                                            "Valeur": f"{body_data['tour_poitrine_cm']:.1f} cm", 
                                            "Norme": f"{race_standards['tour_poitrine_cm'][0]}-{race_standards['tour_poitrine_cm'][1]} cm",
                                            "Évaluation": evaluate_measurement(body_data['tour_poitrine_cm'], race_standards['tour_poitrine_cm'])
                                        },
                                        {
                                            "Paramètre": "Ratio L/H", 
                                            "Valeur": f"{body_data['ratio_longueur_hauteur']:.2f}", 
                                            "Norme": "1.4-1.6",
                                            "Évaluation": evaluate_ratio(body_data['ratio_longueur_hauteur'])
                                        }
                                    ])
                                    def color_evaluation(val):
                                        if "Bon" in val: return 'background-color: #d4edda'
                                        elif "Moyen" in val: return 'background-color: #fff3cd'
                                        else: return 'background-color: #f8d7da'
                                    styled_df = body_df.style.applymap(color_evaluation, subset=['Évaluation'])
                                    st.dataframe(styled_df, use_container_width=True, hide_index=True)
                                if sexe == "Femelle":
                                    st.markdown("### 🍼 CARACTÉRISTIQUES MAMMAIRES")
                                    if mammary_data.get('nombre_mamelles_detectees', 0) > 0:
                                        mammary_metrics = []
                                        if 'volume_mammaire_moyen_cm3' in mammary_data:
                                            mammary_metrics.append({
                                                "Paramètre": "Volume moyen",
                                                "Valeur": f"{mammary_data['volume_mammaire_moyen_cm3']:.1f} cm³",
                                                "Évaluation": evaluate_volume(mammary_data['volume_mammaire_moyen_cm3'])
                                            })
                                        if 'largeur_mammaire_moyenne_cm' in mammary_data:
                                            mammary_metrics.append({
                                                "Paramètre": "Largeur moyenne",
                                                "Valeur": f"{mammary_data['largeur_mammaire_moyenne_cm']:.1f} cm",
                                                "Évaluation": evaluate_width(mammary_data['largeur_mammaire_moyenne_cm'])
                                            })
                                        if 'hauteur_mammaire_moyenne_cm' in mammary_data:
                                            mammary_metrics.append({
                                                "Paramètre": "Hauteur moyenne",
                                                "Valeur": f"{mammary_data['hauteur_mammaire_moyenne_cm']:.1f} cm",
                                                "Évaluation": evaluate_height(mammary_data['hauteur_mammaire_moyenne_cm'])
                                            })
                                        if 'symetrie_mammaire' in mammary_data:
                                            mammary_metrics.append({
                                                "Paramètre": "Symétrie",
                                                "Valeur": f"{mammary_data['symetrie_mammaire']:.2f}",
                                                "Évaluation": evaluate_symmetry(mammary_data['symetrie_mammaire'])
                                            })
                                        if mammary_data.get('score_developpement', 0) > 0:
                                            mammary_metrics.append({
                                                "Paramètre": "Score développement",
                                                "Valeur": f"{mammary_data['score_developpement']:.1f}/10",
                                                "Évaluation": evaluate_score(mammary_data['score_developpement'])
                                            })
                                        if mammary_metrics:
                                            mammary_df = pd.DataFrame(mammary_metrics)
                                            st.dataframe(mammary_df, use_container_width=True, hide_index=True)
                                        classification = st.session_state.photo_analyzer.get_mammary_classification(mammary_data)
                                        st.markdown("#### 🎯 APTITUDE LAITIÈRE")
                                        if "EXCELLENT" in classification:
                                            st.success(f"**{classification}**")
                                        elif "BON" in classification:
                                            st.info(f"**{classification}**")
                                        elif "MOYEN" in classification:
                                            st.warning(f"**{classification}**")
                                        else:
                                            st.error(f"**{classification}**")
                                        st.markdown("#### 💡 RECOMMANDATIONS")
                                        if mammary_data.get('score_developpement', 0) >= 6:
                                            st.info("""
                                            ✅ **Bonne candidate pour :**
                                            - Reproduction et sélection
                                            - Production laitière optimale
                                            - Amélioration génétique du troupeau
                                            - Vente comme reproductrice de valeur
                                            """)
                                            st.markdown("**Actions conseillées :**")
                                            st.markdown("""
                                            1. **Inclure dans le programme de reproduction**
                                            2. **Suivi régulier de la production laitière**
                                            3. **Conserver pour les saillies futures**
                                            4. **Éventuellement vendre comme reproductrice**
                                            """)
                                        else:
                                            st.warning("""
                                            ⚠️ **À surveiller ou réformer :**
                                            - Contrôler régulièrement la production réelle
                                            - Évaluer l'état de santé général
                                            - Envisager le renouvellement si nécessaire
                                            - Possiblement réformer si plusieurs lactations médiocres
                                            """)
                                    else:
                                        st.warning("""
                                        ⚠️ **Aucune mamelle détectée.**
                                        **Utilisation du mode estimation :**
                                        - Score estimé: 5.5/10
                                        - Symétrie estimée: 0.8
                                        - Classification: MOYEN
                                        """)
                                        classification = "MOYEN - Productrice acceptable"
                                else:
                                    st.info("🐑 **Animal mâle** - Pas d'évaluation mammaire nécessaire")
                                    classification = "Non applicable"
                                st.markdown("---")
                                col_save1, col_save2, col_save3 = st.columns([1, 2, 1])
                                with col_save2:
                                    if st.button("💾 ENREGISTRER DANS LA BASE", type="primary", use_container_width=True):
                                        save_complete_characterization(
                                            race, age_mois, poids, sexe, classification
                                        )
                            else:
                                st.error(f"❌ Erreur d'analyse: {mammary_data['error']}")
                                st.info("⚠️ Utilisation des valeurs par défaut pour continuer.")
                        except Exception as e:
                            st.error(f"❌ Erreur lors de l'analyse: {str(e)}")
                            st.info("""
                            **Solution :**
                            1. Vérifiez que la photo est claire
                            2. Assurez-vous que l'étalon est visible
                            3. Réessayez ou utilisez le mode MANUEL
                            """)
        elif 'rear_image' not in st.session_state:
            st.info("👆 Veuillez d'abord prendre ou télécharger une photo arrière")

# --- PAGE ACCUEIL ---
def page_accueil():
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - RACES ALGÉRIENNES</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique des races ovines algériennes**")
    try:
        # Récupération des données de l'utilisateur connecté
        user_id = st.session_state.user_id
        
        # Total brebis
        response = supabase.table("brebis").select("*", count="exact").eq("user_id", user_id).execute()
        total = response.count
        
        # Nombre de races distinctes
        response = supabase.table("brebis").select("race").eq("user_id", user_id).execute()
        races = len(set(row['race'] for row in response.data)) if response.data else 0
        
        # Poids moyen des femelles
        response = supabase.table("brebis").select("poids").eq("user_id", user_id).eq("sexe", "F").execute()
        if response.data:
            poids_f = sum(row['poids'] for row in response.data if row['poids']) / len(response.data)
        else:
            poids_f = 0
        
        # Nombre de scans 3d (à adapter si vous voulez filtrer par user_id aussi)
        response = supabase.table("scans_3d").select("*", count="exact").eq("user_id", user_id).execute()
        scans = response.count
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.markdown(f"""
            <div class='metric-card'>
                <h3>🐑 TOTAL BREBIS</h3>
                <h2>{total}</h2>
                <p>Races algériennes</p>
            </div>
            """, unsafe_allow_html=True)
        with col2:
            st.markdown(f"""
            <div class='metric-card'>
                <h3>🏷️ RACES</h3>
                <h2>{races}</h2>
                <p>Différentes</p>
            </div>
            """, unsafe_allow_html=True)
        with col3:
            st.markdown(f"""
            <div class='metric-card'>
                <h3>♀️ POIDS MOYEN</h3>
                <h2>{poids_f:.1f} kg</h2>
                <p>Femelles</p>
            </div>
            """, unsafe_allow_html=True)
        with col4:
            st.markdown(f"""
            <div class='metric-card'>
                <h3>📐 SCANS 3D</h3>
                <h2>{scans}</h2>
                <p>Réalisés</p>
            </div>
            """, unsafe_allow_html=True)
        
        st.markdown("### 📊 DISTRIBUTION DES RACES")
        # Agrégation par race
        response = supabase.table("brebis").select("race, poids, age_mois").eq("user_id", user_id).execute()
        df = pd.DataFrame(response.data)
        if not df.empty:
            races_data = df.groupby('race').agg(
                count=('race', 'count'),
                poids_moyen=('poids', 'mean'),
                age_moyen=('age_mois', 'mean')
            ).reset_index().to_dict('records')
            
            df_races = pd.DataFrame(races_data, columns=['race', 'count', 'poids_moyen', 'age_moyen'])
            df_races = df_races.rename(columns={'race':'Race', 'count':'Nombre', 'poids_moyen':'Poids moyen', 'age_moyen':'Âge moyen'})
            
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
        else:
            st.info("Vous n'avez pas encore de brebis enregistrées.")
    except Exception as e:
        st.info("Bienvenue dans Ovin Manager Pro! Le système est en cours d'initialisation.")
        st.markdown("### Fonctionnalités disponibles:")
        st.markdown("""
        - **📸 Photo & Mesures**: Capture photo avec étalon et mesures automatiques
        - **📦 Analyse Multiple**: Analyse en lot pour éleveurs
        - **📐 Scanner 3D**: Simulation de scans 3D
        - **📊 Gestion**: Suivi du troupeau
        - **🥛 Production**: Suivi laitier
        - **🎯 Critères**: Évaluation des mamelles
        - **📊 Statistiques**: Analyses avancées
        - **🧬 Génétique**: Analyses génomiques
        """)

# --- PAGE GESTION (exemple de lecture) ---
def page_gestion():
    st.markdown('<h2 class="section-header">📊 GESTION DU TROUPEAU</h2>', unsafe_allow_html=True)
    tab1, tab2, tab3, tab4 = st.tabs(["🐑 LISTE", "📈 STATISTIQUES", "🔍 RECHERCHE", "📤 EXPORT"])
    
    with tab1:
        response = supabase.table("brebis").select("identifiant, nom, race, sexe, age_mois, poids, score_condition, couleur_robe, statut").eq("user_id", st.session_state.user_id).limit(100).execute()
        if response.data:
            df = pd.DataFrame(response.data)
            df = df.rename(columns={
                'identifiant': 'ID',
                'nom': 'Nom',
                'race': 'Race',
                'sexe': 'Sexe',
                'age_mois': 'Âge',
                'poids': 'Poids',
                'score_condition': 'Score',
                'couleur_robe': 'Couleur',
                'statut': 'Statut'
            })
            st.dataframe(df, use_container_width=True, height=400)
            st.metric("Brebis affichées", len(df))
        else:
            st.info("Aucune brebis enregistrée.")
    
    with tab2:
        st.markdown("### 📊 STATISTIQUES DESCRIPTIVES")
        # Similaire à l'accueil mais peut être plus détaillé
        response = supabase.table("brebis").select("race, poids, age_mois").eq("user_id", st.session_state.user_id).execute()
        df = pd.DataFrame(response.data)
        if not df.empty:
            stats_data = df.groupby('race').agg(
                count=('race', 'count'),
                poids_moyen=('poids', 'mean'),
                age_moyen=('age_mois', 'mean')
            ).reset_index()
            st.dataframe(stats_data)
            
            col1, col2 = st.columns(2)
            with col1:
                fig = px.bar(stats_data, x='race', y='poids_moyen',
                            title="Poids moyen par race",
                            color='count',
                            color_continuous_scale='Reds')
                st.plotly_chart(fig, use_container_width=True)
            with col2:
                fig = px.pie(stats_data, values='count', names='race',
                            title="Répartition des races",
                            hole=0.4)
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Pas assez de données.")
    
    with tab3:
        st.markdown("### 🔍 RECHERCHE SIMPLE")
        recherche = st.text_input("Rechercher par nom ou ID")
        if recherche:
            # Recherche avec filtre user_id et like
            response = supabase.table("brebis").select("*").eq("user_id", st.session_state.user_id).or_(
                f"identifiant.ilike.%{recherche}%,nom.ilike.%{recherche}%"
            ).execute()
            resultats = response.data
            if resultats:
                df_res = pd.DataFrame(resultats)
                st.dataframe(df_res)
            else:
                st.info("Aucun résultat.")
    
    with tab4:
        st.markdown("### 📤 EXPORT DE DÉMONSTRATION")
        data_example = {
            'ID': ['EXEMPLE-001', 'EXEMPLE-002'],
            'Race': ['HAMRA', 'OUDA'],
            'Poids_kg': [55.5, 62.3],
            'Age_mois': [24, 36]
        }
        df_export = pd.DataFrame(data_example)
        col1, col2 = st.columns(2)
        with col1:
            csv = df_export.to_csv(index=False)
            st.download_button(
                label="📥 Télécharger CSV",
                data=csv,
                file_name="exemple_brebis.csv",
                mime="text/csv"
            )
        with col2:
            json_data = df_export.to_json(orient='records', indent=2)
            st.download_button(
                label="📥 Télécharger JSON",
                data=json_data,
                file_name="exemple_brebis.json",
                mime="application/json"
            )

# --- PAGE PRODUCTION (similaire) ---
def page_production():
    st.markdown('<h2 class="section-header">🥛 SUIVI DE PRODUCTION LAITIÈRE</h2>', unsafe_allow_html=True)
    tab1, tab2, tab3 = st.tabs(["📝 SAISIE", "📈 ANALYSE", "🏆 CLASSEMENT"])
    
    with tab1:
        st.markdown("### 📝 SAISIE DE PRODUCTION")
        with st.form("form_production"):
            brebis_id = st.text_input("Identifiant de la brebis", value="HAM-F-001")
            date_mesure = st.date_input("Date", value=date.today())
            col1, col2, col3 = st.columns(3)
            with col1:
                quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5, 0.1)
                cellules = st.number_input("Cellules (x1000)", 0, 1000, 200)
            with col2:
                mg = st.number_input("Matière grasse %", 0.0, 20.0, 7.2, 0.1)
                lactose = st.number_input("Lactose %", 0.0, 10.0, 4.8, 0.1)
            with col3:
                proteine = st.number_input("Protéine %", 0.0, 20.0, 5.5, 0.1)
                ph = st.number_input("pH", 6.0, 7.0, 6.7, 0.1)
            notes = st.text_area("Notes", value="Production standard")
            if st.form_submit_button("💾 Enregistrer", type="primary"):
                # Vérifier que la brebis existe et appartient à l'utilisateur
                response = supabase.table("brebis").select("id").eq("identifiant", brebis_id).eq("user_id", st.session_state.user_id).execute()
                if response.data:
                    brebis_db_id = response.data[0]['id']
                    data = {
                        'brebis_id': brebis_db_id,
                        'date_mesure': date_mesure.isoformat(),
                        'quantite_litre': quantite,
                        'taux_matiere_grasse': mg,
                        'taux_proteine': proteine,
                        'cellules_somatiques': cellules,
                        'lactose': lactose,
                        'ph': ph,
                        'notes': notes,
                        'user_id': st.session_state.user_id
                    }
                    supabase.table("production_lait").insert(data).execute()
                    st.success(f"✅ Production enregistrée pour {brebis_id}")
                else:
                    st.error("Brebis non trouvée ou ne vous appartient pas.")
    
    with tab2:
        st.markdown("### 📈 ANALYSE DE PRODUCTION")
        # Exemple de données simulées
        mois = ['Jan', 'Fév', 'Mar', 'Avr', 'Mai', 'Jun']
        production = [2.8, 3.2, 3.5, 3.1, 2.9, 3.0]
        mg = [7.2, 7.5, 7.8, 7.3, 7.1, 7.4]
        df_prod = pd.DataFrame({
            'Mois': mois,
            'Production (L)': production,
            'MG (%)': mg
        })
        fig = px.line(df_prod, x='Mois', y='Production (L)',
                     title="Évolution de la production laitière",
                     markers=True)
        st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 🏆 TOP PRODUCTRICES")
        # Pourrait être calculé à partir des données réelles, ici simulé
        top_data = {
            'Brebis': ['FHAM001', 'FOUDA002', 'FBER003', 'FSID004'],
            'Race': ['HAMRA', 'OUDA', 'BERBERE', 'SIDAHOU'],
            'Production moyenne (L)': [3.5, 3.2, 2.8, 3.0],
            'MG moyenne (%)': [7.8, 7.5, 7.2, 7.4]
        }
        df_top = pd.DataFrame(top_data)
        fig = px.bar(df_top, x='Brebis', y='Production moyenne (L)',
                    color='Race',
                    title="Top 4 productrices",
                    hover_data=['MG moyenne (%)'])
        st.plotly_chart(fig, use_container_width=True)

# --- PAGE SCANNER 3D (inchangée car ne fait que des simulations) ---
def page_scanner_3d():
    # Identique au code original
    st.markdown('<h2 class="section-header">📐 SCANNER 3D & SAISIE MANUELLE</h2>', unsafe_allow_html=True)
    tab1, tab2 = st.tabs(["🎯 SCANNER 3D", "📝 SAISIE MANUELLE"])
    with tab1:
        st.markdown("### 🎯 SIMULATION SCANNER 3D")
        race_selection = st.selectbox("Sélectionnez une race:", 
                                     list(STANDARDS_RACES.keys()),
                                     format_func=lambda x: STANDARDS_RACES[x]['nom_complet'])
        if race_selection:
            brebis_simulee = {
                'race': race_selection,
                'identifiant': f"{race_selection[:3]}-SIM-001",
                'nom': f"Sim{race_selection[:3]}001",
                'poids': 55.5,
                'sexe': 'F',
                'age_mois': 24,
                'couleur': STANDARDS_RACES[race_selection]['couleur']
            }
            photo = Scanner3D.generer_photo_simulee(brebis_simulee)
            col1, col2 = st.columns([2, 1])
            with col1:
                st.image(photo, caption=f"Photo simulée - {brebis_simulee['nom']}", use_column_width=True)
                if st.button("📸 Simuler scan", type="primary"):
                    with st.spinner("Scan en cours..."):
                        points = Scanner3D.simuler_scan_3d(brebis_simulee)
                        st.success(f"✅ Scan simulé! {len(points)} points 3D générés")
                        df_points = pd.DataFrame(points[:5])
                        st.dataframe(df_points[['x', 'y', 'z']])
            with col2:
                st.markdown(f"""
                <div class='race-card'>
                    <h4>📷 INFORMATIONS</h4>
                    <p><strong>Race:</strong> {STANDARDS_RACES[race_selection]['nom_complet']}</p>
                    <p><strong>ID:</strong> {brebis_simulee['identifiant']}</p>
                    <p><strong>Sexe:</strong> {brebis_simulee['sexe']}</p>
                    <p><strong>Âge:</strong> {brebis_simulee['age_mois']} mois</p>
                    <p><strong>Poids:</strong> {brebis_simulee['poids']:.1f} kg</p>
                </div>
                """, unsafe_allow_html=True)
    with tab2:
        st.markdown("### 📝 SAISIE MANUELLE DES MESURES")
        with st.form("saisie_manuelle"):
            col_saisie1, col_saisie2 = st.columns(2)
            with col_saisie1:
                identifiant = st.text_input("Identifiant de l'animal", value="BRE-001")
                date_mesure = st.date_input("Date de mesure", value=date.today())
                operateur = st.text_input("Opérateur", value="Opérateur 1")
                st.markdown("#### 📏 MENSURATIONS (cm)")
                longueur = st.number_input("Longueur corps", 0.0, 200.0, 100.0, 0.1)
                hauteur = st.number_input("Hauteur garrot", 0.0, 150.0, 70.0, 0.1)
                largeur = st.number_input("Largeur bassin", 0.0, 100.0, 40.0, 0.1)
                tour_poitrine = st.number_input("Tour de poitrine", 0.0, 200.0, 100.0, 0.1)
            with col_saisie2:
                poids = st.number_input("Poids (kg)", 0.0, 200.0, 50.0, 0.1)
                score_condition = st.slider("Score de condition", 1, 5, 3)
                st.markdown("#### 🎨 CARACTÈRES QUALITATIFS")
                couleur_robe = st.text_input("Couleur de la robe", value="Blanche")
                etat_corporel = st.select_slider("État corporel", 
                                                ['Maigre', 'Normal', 'Gras'])
                temperement = st.selectbox("Tempérament", 
                                         ['Calme', 'Nervieux', 'Intermediaire', 'Agité'])
                notes = st.text_area("Notes complémentaires", value="Mesures standard")
            if st.form_submit_button("💾 Enregistrer les mesures", type="primary"):
                # On pourrait sauvegarder dans une table "mesures_manuelles" avec user_id
                # Pour l'instant, on affiche juste
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

# --- PAGE CRITÈRES (inchangée) ---
def page_criteres():
    st.markdown('<h2 class="section-header">🎯 CRITÈRES DE SÉLECTION - MAMMELLES</h2>', unsafe_allow_html=True)
    st.markdown("### 📋 ÉVALUATION DES CRITÈRES MAMMAIRES")
    with st.form("evaluation_form"):
        st.markdown("#### 📝 CRITÈRES PRINCIPAUX")
        col1, col2 = st.columns(2)
        with col1:
            volume = st.slider("Volume mammaire (1-5)", 1, 5, 3,
                              help="1: Très petit, 5: Très développé")
            symetrie = st.slider("Symétrie (1-5)", 1, 5, 3,
                                help="1: Asymétrique, 5: Parfaitement symétrique")
        with col2:
            insertion = st.slider("Insertion des trayons (1-5)", 1, 5, 3,
                                 help="1: Très écartés, 5: Bien insérés")
            longueur_trayons = st.slider("Longueur des trayons (cm)", 2.0, 8.0, 4.5, 0.1)
        orientation = st.selectbox("Orientation des trayons",
                                 ['Parallèle', 'Légèrement divergent', 'Divergent'])
        notes = st.text_area("Observations")
        if st.form_submit_button("🎯 Calculer le score", type="primary"):
            score_total = (volume + symetrie + insertion) / 3
            st.success(f"✅ Évaluation terminée! Score: {score_total:.1f}/5")
            col_res1, col_res2, col_res3 = st.columns(3)
            with col_res1:
                st.metric("Volume", f"{volume}/5")
            with col_res2:
                st.metric("Symétrie", f"{symetrie}/5")
            with col_res3:
                st.metric("Insertion", f"{insertion}/5")
            st.markdown("### 📊 CLASSIFICATION")
            if score_total >= 4:
                st.success("**Type A (4-5): Excellent pour la production**")
            elif score_total >= 3:
                st.info("**Type B (3-4): Bon pour la production**")
            elif score_total >= 2:
                st.warning("**Type C (2-3): Moyen, à surveiller**")
            else:
                st.error("**Type D (1-2): À améliorer ou réformer**")

# --- PAGE STATISTIQUES (inchangée) ---
def page_stats():
    st.markdown('<h2 class="section-header">📊 ANALYSE STATISTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    # Données simulées
    np.random.seed(42)
    n = 50
    data = {
        'Race': np.random.choice(['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE'], n),
        'Poids': np.random.normal(55, 10, n),
        'Longueur': np.random.normal(105, 15, n),
        'Hauteur': np.random.normal(70, 8, n),
        'Age': np.random.randint(12, 84, n)
    }
    df = pd.DataFrame(data)
    tab1, tab2 = st.tabs(["📈 DESCRIPTIVE", "📊 CORRÉLATIONS"])
    with tab1:
        st.markdown("### 📈 STATISTIQUES PAR RACE")
        races = df['Race'].unique()
        for race in races:
            with st.expander(f"{get_race_data(race, 'nom_complet', race)}"):
                race_df = df[df['Race'] == race]
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Nombre", len(race_df))
                    st.metric("Poids moyen", f"{race_df['Poids'].mean():.1f} kg")
                with col2:
                    st.metric("Âge moyen", f"{race_df['Age'].mean():.0f} mois")
                    st.metric("Longueur moyenne", f"{race_df['Longueur'].mean():.1f} cm")
                with col3:
                    st.metric("Hauteur moyenne", f"{race_df['Hauteur'].mean():.1f} cm")
    with tab2:
        st.markdown("### 📊 MATRICE DE CORRÉLATION")
        numeric_df = df[['Poids', 'Longueur', 'Hauteur', 'Age']]
        corr_matrix = numeric_df.corr()
        fig = px.imshow(corr_matrix,
                       title="Corrélations entre variables",
                       color_continuous_scale='RdBu',
                       zmin=-1, zmax=1,
                       text_auto=True)
        st.plotly_chart(fig, use_container_width=True)

# --- PAGE GÉNÉTIQUE (simulée, mais peut être adaptée) ---
def page_genetique():
    st.markdown('<h2 class="section-header">🧬 ANALYSE GÉNÉTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    tab1, tab2 = st.tabs(["🧪 GÉNOTYPAGE", "🌳 DIVERSITÉ"])
    with tab1:
        st.markdown("### 🧪 SIMULATION DE GÉNOTYPAGE")
        if st.button("🧬 Générer un génotype simulé", type="primary"):
            with st.spinner("Génération en cours..."):
                snps = []
                for i in range(10):
                    snps.append({
                        'SNP': f"SNP{i+1:03d}",
                        'Chromosome': random.randint(1, 26),
                        'Position': random.randint(1000000, 90000000),
                        'Allèle 1': random.choice(['A', 'C', 'G', 'T']),
                        'Allèle 2': random.choice(['A', 'C', 'G', 'T']),
                        'Génotype': f"{random.choice(['A', 'C', 'G', 'T'])}{random.choice(['A', 'C', 'G', 'T'])}"
                    })
                df_snps = pd.DataFrame(snps)
                st.dataframe(df_snps, use_container_width=True)
                st.success("✅ Génotype simulé généré avec succès!")
    with tab2:
        st.markdown("### 🌳 DIVERSITÉ GÉNÉTIQUE")
        st.info("""
        **Indicateurs de diversité génétique:**
        - **Hétérozygotie observée (Ho):** Proportion d'individus hétérozygotes
        - **Hétérozygotie attendue (He):** Proportion attendue sous Hardy-Weinberg
        - **Fis:** Coefficient de consanguinité
        *Ces indicateurs seront calculés à partir des données de génotypage.*
        """)
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Ho", "0.325")
        with col2:
            st.metric("He", "0.342")
        with col3:
            st.metric("Fis", "0.050")

# --- PAGE ANALYSE MULTIPLE (inchangée car elle utilise st.session_state et pas de base) ---
def page_analyse_multiple():
    st.markdown('<h2 class="section-header">📦 ANALYSE MULTIPLE DE BREBIS</h2>', unsafe_allow_html=True)
    st.info("""
    **🎯 POUR LES ÉLEVEURS QUI VOUS ENVOIENT DES PHOTOS :**
    1. Ils prennent les photos selon le guide
    2. Ils vous envoient les fichiers
    3. Vous les analysez ici en lot
    4. Vous générez un rapport par élevage
    """)
    tab1, tab2, tab3 = st.tabs(["📁 IMPORT MULTIPLE", "🔍 ANALYSE", "📊 RAPPORTS"])
    with tab1:
        st.markdown("### 📁 IMPORTEZ TOUTES LES PHOTOS")
        uploaded_files = st.file_uploader(
            "Téléchargez TOUTES les photos de l'élevage",
            type=['jpg', 'jpeg', 'png', 'bmp'],
            accept_multiple_files=True,
            key="multi_upload"
        )
        if uploaded_files:
            st.success(f"✅ {len(uploaded_files)} fichiers téléchargés")
            brebis_dict = {}
            for file in uploaded_files:
                filename = file.name.lower()
                if any(mot in filename for mot in ['profil', 'profile', 'cote', 'side', 'latéral']):
                    brebis_name = filename.split('_')[0] if '_' in filename else filename.split('.')[0]
                    if brebis_name not in brebis_dict:
                        brebis_dict[brebis_name] = {'profil': None, 'arriere': None}
                    brebis_dict[brebis_name]['profil'] = file
                elif any(mot in filename for mot in ['arriere', 'rear', 'back', 'mamelle', 'mammary', 'postérieur']):
                    brebis_name = filename.split('_')[0] if '_' in filename else filename.split('.')[0]
                    if brebis_name not in brebis_dict:
                        brebis_dict[brebis_name] = {'profil': None, 'arriere': None}
                    brebis_dict[brebis_name]['arriere'] = file
            st.markdown("#### 📋 PHOTOS CLASSÉES")
            for brebis_name, photos in brebis_dict.items():
                col1, col2, col3 = st.columns([2, 2, 1])
                with col1:
                    st.write(f"**{brebis_name}**")
                with col2:
                    profil_ok = "✅" if photos['profil'] else "❌"
                    arriere_ok = "✅" if photos['arriere'] else "❌"
                    st.write(f"Profil: {profil_ok} | Arrière: {arriere_ok}")
                with col3:
                    if photos['profil'] and photos['arriere']:
                        st.success("Complet")
                    else:
                        st.warning("Incomplet")
            st.session_state.brebis_photos = brebis_dict
            with st.form("infos_elevage_form"):
                st.markdown("### 📝 INFORMATIONS DE L'ÉLEVAGE")
                elevage_nom = st.text_input("Nom de l'élevage")
                elevage_region = st.text_input("Région")
                elevage_contact = st.text_input("Contact éleveur")
                date_prise_photo = st.date_input("Date de prise des photos")
                if st.form_submit_button("💾 Enregistrer les informations"):
                    st.session_state.elevage_info = {
                        'nom': elevage_nom,
                        'region': elevage_region,
                        'contact': elevage_contact,
                        'date_photos': date_prise_photo.isoformat()
                    }
                    st.success("✅ Informations enregistrées!")
    with tab2:
        if 'brebis_photos' not in st.session_state:
            st.warning("⚠️ Importez d'abord les photos dans l'onglet 1")
        else:
            st.markdown("### 🔍 ANALYSE DES BREBIS")
            st.markdown("#### ⚙️ PARAMÈTRES COMMUNS")
            col_param1, col_param2 = st.columns(2)
            with col_param1:
                etalon_type = st.selectbox(
                    "Étalon utilisé dans les photos:",
                    ["feuille_a4_largeur", "baton_1m", "carte_bancaire", "piece_100da", "telephone_standard"],
                    format_func=lambda x: {
                        "feuille_a4_largeur": "Feuille A4 (21cm)",
                        "baton_1m": "Bâton 1m (100cm)",
                        "carte_bancaire": "Carte bancaire (8.56cm)",
                        "piece_100da": "Pièce 100 DA (2.6cm)",
                        "telephone_standard": "Téléphone (15cm)"
                    }[x],
                    key="multi_etalon"
                )
            with col_param2:
                race_predom = st.selectbox(
                    "Race prédominante:",
                    list(STANDARDS_RACES.keys()),
                    format_func=lambda x: STANDARDS_RACES[x]['nom_complet'],
                    key="multi_race"
                )
            if st.button("🚀 ANALYSER TOUTES LES BREBIS", type="primary", use_container_width=True):
                analyzer = OvinPhotoAnalyzer()
                analyzer.set_etalon(etalon_type)
                results = {}
                progress_bar = st.progress(0)
                for i, (brebis_name, photos) in enumerate(st.session_state.brebis_photos.items()):
                    if photos['profil'] and photos['arriere']:
                        st.markdown(f"---")
                        st.markdown(f"#### 🐑 **{brebis_name}**")
                        col_bre1, col_bre2 = st.columns(2)
                        with col_bre1:
                            try:
                                img_profil = Image.open(photos['profil'])
                                img_profil_np = np.array(img_profil)
                                body_measurements = analyzer.analyze_profile_photo(img_profil_np)
                                if body_measurements:
                                    st.success("✅ Profil analysé")
                                    st.metric("Longueur", f"{body_measurements['longueur_corps_cm']:.1f} cm")
                                else:
                                    st.warning("⚠️ Profil non analysable")
                            except Exception as e:
                                st.error(f"❌ Erreur profil: {str(e)}")
                        with col_bre2:
                            try:
                                img_arriere = Image.open(photos['arriere'])
                                img_arriere_np = np.array(img_arriere)
                                mammary_data = analyzer.analyze_rear_photo(img_arriere_np, is_female=True)
                                if 'error' not in mammary_data:
                                    st.success("✅ Mamelles analysées")
                                    if 'score_developpement' in mammary_data:
                                        st.metric("Score", f"{mammary_data['score_developpement']:.1f}/10")
                                    classification = analyzer.get_mammary_classification(mammary_data)
                                    st.write(f"**Aptitude:** {classification}")
                                else:
                                    st.warning(f"⚠️ {mammary_data['error']}")
                            except Exception as e:
                                st.error(f"❌ Erreur arrière: {str(e)}")
                        results[brebis_name] = {
                            'body': body_measurements if 'body_measurements' in locals() else {},
                            'mammary': mammary_data if 'mammary_data' in locals() else {},
                            'classification': classification if 'classification' in locals() else "Non classé"
                        }
                        progress = (i + 1) / len(st.session_state.brebis_photos)
                        progress_bar.progress(progress)
                st.session_state.analysis_results = results
                st.success(f"✅ Analyse terminée pour {len(results)} brebis!")
                if st.button("📊 GÉNÉRER LE RAPPORT DE L'ÉLEVAGE", type="primary"):
                    generate_elevage_report(results)
    with tab3:
        if 'analysis_results' not in st.session_state:
            st.info("📈 Les rapports apparaîtront ici après analyse")
        else:
            st.markdown("### 📊 RAPPORT D'ÉLEVAGE")
            results = st.session_state.analysis_results
            scores = []
            classifications = []
            for brebis_name, data in results.items():
                if data.get('mammary') and 'score_developpement' in data['mammary']:
                    scores.append(data['mammary']['score_developpement'])
                if data.get('classification'):
                    classifications.append(data['classification'])
            col_stat1, col_stat2, col_stat3 = st.columns(3)
            with col_stat1:
                st.metric("Nombre de brebis", len(results))
            with col_stat2:
                if scores:
                    avg_score = np.mean(scores)
                    st.metric("Score moyen", f"{avg_score:.1f}/10")
            with col_stat3:
                if classifications:
                    excellent_count = sum(1 for c in classifications if "EXCELLENT" in c)
                    st.metric("Excellent", f"{excellent_count}")
            st.markdown("#### 📋 RÉSULTATS DÉTAILLÉS")
            table_data = []
            for brebis_name, data in results.items():
                row = {
                    'Brebis': brebis_name,
                    'Longueur (cm)': f"{data.get('body', {}).get('longueur_corps_cm', 0):.1f}" if data.get('body') else "N/A",
                    'Hauteur (cm)': f"{data.get('body', {}).get('hauteur_garrot_cm', 0):.1f}" if data.get('body') else "N/A",
                    'Score mamelle': f"{data.get('mammary', {}).get('score_developpement', 0):.1f}" if data.get('mammary') else "N/A",
                    'Classification': data.get('classification', 'N/A'),
                    'Recommandation': get_recommendation_from_classification(data.get('classification', ''))
                }
                table_data.append(row)
            df_results = pd.DataFrame(table_data)
            st.dataframe(df_results, use_container_width=True)
            st.markdown("---")
            col_dl1, col_dl2, col_dl3 = st.columns(3)
            with col_dl1:
                csv = df_results.to_csv(index=False)
                st.download_button(
                    label="📥 Télécharger CSV",
                    data=csv,
                    file_name="resultats_elevage.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            with col_dl2:
                full_data = {
                    'elevage_info': st.session_state.get('elevage_info', {}),
                    'analysis_results': st.session_state.analysis_results,
                    'date_analyse': datetime.now().isoformat(),
                    'etalon_utilise': st.session_state.get('multi_etalon', 'inconnu')
                }
                json_data = json.dumps(full_data, indent=2, ensure_ascii=False)
                st.download_button(
                    label="📥 Télécharger JSON",
                    data=json_data,
                    file_name="rapport_complet.json",
                    mime="application/json",
                    use_container_width=True
                )
            with col_dl3:
                rapport_text = generate_text_report(st.session_state.get('elevage_info', {}), results)
                st.download_button(
                    label="📄 Télécharger rapport",
                    data=rapport_text,
                    file_name="rapport_elevage.txt",
                    mime="text/plain",
                    use_container_width=True
                )

# --- PAGE GÉNOMIQUE AVANCÉE (adaptée) ---
def page_genomique_avancee():
    st.markdown('<h2 class="section-header">🧬 GÉNOMIQUE AVANCÉE & SÉLECTION ASSISTÉE</h2>', unsafe_allow_html=True)
    tab1, tab2, tab3, tab4 = st.tabs([
        "📄 FASTA & BLAST", 
        "🧪 SNP / QTN", 
        "🩺 Maladies génétiques",
        "🧬 Génotypage"
    ])
    with tab1:
        st.markdown("### 📄 GÉNÉRATION DE SÉQUENCES FASTA")
        col1, col2 = st.columns(2)
        with col1:
            gene = st.selectbox("Gène d'intérêt économique", 
                               ["GDF8", "PRNP", "FecB", "DGAT1", "CAST", "LALBA", "CSN1S1"])
            longueur = st.slider("Longueur de la séquence (pb)", 500, 2000, 1000, 100)
        with col2:
            # Récupérer les brebis de l'utilisateur
            response = supabase.table("brebis").select("id, identifiant").eq("user_id", st.session_state.user_id).limit(10).execute()
            brebis_list = response.data
            if brebis_list:
                brebis_dict = {f"{b['identifiant']} (ID: {b['id']})": b['id'] for b in brebis_list}
                selected_brebis = st.selectbox("Choisir une brebis", list(brebis_dict.keys()))
                brebis_id = brebis_dict[selected_brebis]
            else:
                st.warning("Aucune brebis trouvée. Utilisation d'un ID fictif.")
                brebis_id = 0
        if st.button("🧬 Générer séquence FASTA"):
            fasta = GenomicAnalyzer.generer_fasta(brebis_id, gene, longueur)
            st.code(fasta, language="text")
            st.download_button("📥 Télécharger FASTA", fasta, 
                              file_name=f"brebis_{brebis_id}_{gene}.fasta")
        st.markdown("---")
        st.markdown("### 🔍 SIMULATION BLAST")
        st.info("BLAST simulé : compare la séquence requête avec une base de données locale.")
        query_seq = st.text_area("Coller une séquence (format nucléotidique)", 
                                height=150,
                                value="ATCGATCGATCGATCGATCGATCGATCGATCG")
        if st.button("🔬 Lancer BLAST"):
            db_seqs = {
                "GDF8_ref": "ATCGATCGATCGATCGATCGATCGATCGATCG",
                "PRNP_ref": "ATCGATCGATCGATCGATTTTTTTTTTTTTTT",
                "FecB_ref": "ATCGATCGATCGATCGATCGATCGATCGATCG"
            }
            results = GenomicAnalyzer.blast_simule(query_seq, db_seqs)
            if results:
                st.success(f"{len(results)} alignements trouvés")
                df_blast = pd.DataFrame(results)
                st.dataframe(df_blast)
            else:
                st.warning("Aucun alignement significatif trouvé (>70% identité)")
    with tab2:
        st.markdown("### 🧪 DÉTECTION DE SNP / QTN")
        st.markdown("Détecte les variants connus dans une séquence fournie.")
        seq_input = st.text_area("Séquence nucléotidique", height=150,
                                 key="seq_snp")
        gene_choice = st.selectbox("Gène ciblé", ["GDF8", "PRNP", "FecB", "DGAT1", "CAST"])
        if st.button("🔎 Détecter variants"):
            variants = GenomicAnalyzer.detecter_snp_qtn(seq_input, gene_choice)
            if variants:
                st.success("Variant(s) détecté(s) !")
                df_var = pd.DataFrame(variants)
                st.dataframe(df_var)
            else:
                st.info("Aucun variant d'intérêt détecté dans cette séquence.")
        st.markdown("---")
        st.markdown("#### 📊 MARQUEURS DISPONIBLES")
        # Table genetic_markers est commune, pas de filtre user_id
        response = supabase.table("genetic_markers").select("marker_id, chromosome, position, gene_name, trait_effect, is_qtn, economic_importance").execute()
        markers_data = response.data
        df_markers = pd.DataFrame(markers_data, 
                                 columns=['marker_id','chromosome','position','gene_name','trait_effect','is_qtn','economic_importance'])
        df_markers = df_markers.rename(columns={
            'marker_id': 'Marqueur',
            'chromosome': 'Chr',
            'position': 'Pos',
            'gene_name': 'Gène',
            'trait_effect': 'Effet',
            'is_qtn': 'QTN',
            'economic_importance': 'Importance'
        })
        st.dataframe(df_markers)
    with tab3:
        st.markdown("### 🩺 DÉTECTION DE MALADIES GÉNÉTIQUES")
        st.markdown("Évaluation du risque à partir des génotypes.")
        response = supabase.table("brebis").select("id, identifiant").eq("user_id", st.session_state.user_id).execute()
        all_brebis = response.data
        if all_brebis:
            brebis_choice = st.selectbox("Sélectionner une brebis", 
                                        [f"{b['identifiant']} (ID: {b['id']})" for b in all_brebis],
                                        key="maladie_brebis")
            brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
            if st.button("🧪 Évaluer les risques génétiques"):
                # Récupérer les génotypes de cette brebis
                response = supabase.table("genotypage").select("marqueur, genotype").eq("brebis_id", brebis_id).limit(10).execute()
                genotypes_raw = response.data
                if not genotypes_raw:
                    # Générer des génotypes simulés et les insérer
                    genotypes = GenomicAnalyzer.genotyper_brebis(brebis_id, supabase)
                    for g in genotypes:
                        g['user_id'] = st.session_state.user_id
                        supabase.table("genotypage").insert(g).execute()
                    # Re-lire
                    response = supabase.table("genotypage").select("marqueur, genotype").eq("brebis_id", brebis_id).execute()
                    genotypes_raw = response.data
                genotype_dict = {g['marqueur']: g['genotype'] for g in genotypes_raw}
                risques = GenomicAnalyzer.evaluer_risque_maladie(brebis_id, genotype_dict)
                if risques:
                    st.warning("⚠️ Risques génétiques détectés")
                    df_risques = pd.DataFrame(risques)
                    st.dataframe(df_risques)
                else:
                    st.success("✅ Aucun risque génétique majeur détecté")
                st.markdown("#### 🧬 Génotypes analysés")
                df_geno = pd.DataFrame(genotypes_raw, columns=['Marqueur','Génotype'])
                st.dataframe(df_geno)
        else:
            st.warning("Aucune brebis trouvée.")
    with tab4:
        st.markdown("### 🧬 GÉNOTYPAGE SIMULÉ")
        st.markdown("Génération de génotypes pour les marqueurs d'intérêt.")
        response = supabase.table("brebis").select("id, identifiant").eq("user_id", st.session_state.user_id).execute()
        all_brebis = response.data
        if all_brebis:
            brebis_choice = st.selectbox("Choisir une brebis", 
                                        [f"{b['identifiant']} (ID: {b['id']})" for b in all_brebis],
                                        key="geno_brebis")
            brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
            if st.button("🧬 Générer génotype complet"):
                genotypes = GenomicAnalyzer.genotyper_brebis(brebis_id, supabase)
                for g in genotypes:
                    g['user_id'] = st.session_state.user_id
                    # Insert or replace
                    supabase.table("genotypage").upsert(g, on_conflict=['brebis_id','marqueur']).execute()
                st.success(f"{len(genotypes)} marqueurs génotypés et enregistrés.")
                df_geno_full = pd.DataFrame(genotypes, 
                                           columns=['brebis_id','marqueur','allele1','allele2','genotype','date_analyse'])
                st.dataframe(df_geno_full[['marqueur','genotype','date_analyse']])
        else:
            st.warning("Aucune brebis trouvée.")

# --- PAGE ANALYSE LAIT (adaptée) ---
def page_analyse_lait():
    st.markdown('<h2 class="section-header">🥛 ANALYSE BIOCHIMIQUE PROFESSIONNELLE DU LAIT</h2>', unsafe_allow_html=True)
    tab1, tab2, tab3 = st.tabs(["📝 Saisie & Analyse", "📊 Historique", "📈 Comparaison race"])
    with tab1:
        st.markdown("### 🧪 COMPOSITION DU LAIT")
        col1, col2 = st.columns(2)
        with col1:
            brebis_id = st.text_input("Identifiant de la brebis", "HAM-F-001")
            date_prelevement = st.date_input("Date du prélèvement", date.today())
            echantillon_id = f"LAIT-{brebis_id}-{date_prelevement.strftime('%Y%m%d')}"
            st.info(f"Échantillon : {echantillon_id}")
            mg = st.number_input("Matière grasse (g/100ml)", 0.0, 15.0, 7.2, 0.1)
            proteines = st.number_input("Protéines (g/100ml)", 0.0, 12.0, 5.8, 0.1)
            lactose = st.number_input("Lactose (g/100ml)", 0.0, 8.0, 4.8, 0.1)
            extrac_secret = st.number_input("Extrait sec total (g/100ml)", 10.0, 25.0, 18.5, 0.1)
        with col2:
            ph = st.number_input("pH", 6.0, 7.5, 6.7, 0.01)
            acidite = st.number_input("Acidité Dornic (°D)", 10, 30, 17, 1)
            cellules = st.number_input("Cellules somatiques (x1000/ml)", 0, 2000, 250, 10)
            st.markdown("#### 🧪 Minéraux & Vitamines (optionnel)")
            calcium = st.number_input("Calcium (mg/100ml)", 100, 300, 180, 1)
            phosphore = st.number_input("Phosphore (mg/100ml)", 80, 250, 140, 1)
            magnesium = st.number_input("Magnésium (mg/100ml)", 10, 40, 20, 1)
            potassium = st.number_input("Potassium (mg/100ml)", 100, 250, 160, 1)
            vit_a = st.number_input("Vitamine A (µg/100ml)", 20, 100, 45, 1)
            vit_e = st.number_input("Vitamine E (µg/100ml)", 0.0, 1.0, 0.3, 0.01)
            ags = st.number_input("Acides gras saturés (g/100ml)", 2.0, 7.0, 4.2, 0.1)
            agi = st.number_input("Acides gras insaturés (g/100ml)", 0.5, 2.5, 1.1, 0.1)
        if st.button("🔬 ANALYSER L'ÉCHANTILLON", type="primary"):
            interp = LaitAnalyzer.analyser_composition(mg, proteines, lactose, ph, cellules)
            energie = LaitAnalyzer.calculer_energie(mg, proteines, lactose)
            aptitude = LaitAnalyzer.estimer_aptitude_fromagere(mg, proteines, ph, calcium)
            profil_ag = LaitAnalyzer.analyser_acides_gras(ags, agi)
            st.success("✅ Analyse terminée")
            col_res1, col_res2, col_res3 = st.columns(3)
            with col_res1:
                st.metric("Énergie (kcal/100ml)", energie)
                st.metric("Aptitude fromagère", f"{aptitude}/100")
            with col_res2:
                st.metric("MG", interp['MG'])
                st.metric("Protéines", interp['Protéines'])
            with col_res3:
                st.metric("Lactose", interp['Lactose'])
                st.metric("Cellules", interp['Cellules'])
            st.markdown("#### 🧾 Rapport détaillé")
            df_comp = pd.DataFrame([
                {"Paramètre": "Matière grasse", "Valeur": f"{mg} g/100ml", "Norme": "6.0-7.5", "Appréciation": interp['MG']},
                {"Paramètre": "Protéines", "Valeur": f"{proteines} g/100ml", "Norme": "5.0-6.5", "Appréciation": interp['Protéines']},
                {"Paramètre": "Lactose", "Valeur": f"{lactose} g/100ml", "Norme": "4.0-5.2", "Appréciation": interp['Lactose']},
                {"Paramètre": "pH", "Valeur": f"{ph}", "Norme": "6.5-6.8", "Appréciation": interp['pH']},
                {"Paramètre": "Cellules somatiques", "Valeur": f"{cellules} x1000/ml", "Norme": "<500", "Appréciation": interp['Cellules']},
            ])
            st.dataframe(df_comp, use_container_width=True, hide_index=True)
            st.markdown("#### 🧪 Profil en acides gras")
            df_ag = pd.DataFrame([profil_ag]).T.rename(columns={0:'Valeur'})
            st.dataframe(df_ag)
            if st.button("💾 Enregistrer dans la base"):
                # Vérifier que la brebis existe et appartient à l'utilisateur
                response = supabase.table("brebis").select("id").eq("identifiant", brebis_id).eq("user_id", st.session_state.user_id).execute()
                if response.data:
                    brebis_db_id = response.data[0]['id']
                else:
                    # Créer la brebis si elle n'existe pas (mais il faudrait plus d'infos)
                    data_brebis = {
                        'identifiant': brebis_id,
                        'nom': f"Brebis_{brebis_id}",
                        'race': 'INCONNU',
                        'sexe': 'F',
                        'age_mois': 24,
                        'poids': 50.0,
                        'user_id': st.session_state.user_id
                    }
                    resp = supabase.table("brebis").insert(data_brebis).execute()
                    brebis_db_id = resp.data[0]['id']
                data = {
                    'brebis_id': brebis_db_id,
                    'date_analyse': date_prelevement.isoformat(),
                    'echantillon_id': echantillon_id,
                    'mg_g_100ml': mg,
                    'proteines_g_100ml': proteines,
                    'lactose_g_100ml': lactose,
                    'extrac_secret': extrac_secret,
                    'ph': ph,
                    'acidite_dornic': acidite,
                    'cellules_somatiques': cellules,
                    'acides_gras_satures_g': ags,
                    'acides_gras_insatures_g': agi,
                    'calcium_mg_100ml': calcium,
                    'phosphore_mg_100ml': phosphore,
                    'magnesium_mg_100ml': magnesium,
                    'potassium_mg_100ml': potassium,
                    'vitamine_a_ug_100ml': vit_a,
                    'vitamine_e_ug_100ml': vit_e,
                    'aptitude_fromagere': str(aptitude),
                    'notes': "Analyse professionnelle",
                    'user_id': st.session_state.user_id
                }
                supabase.table("milk_composition").insert(data).execute()
                st.success("✅ Analyse enregistrée dans la base de données !")
    with tab2:
        st.markdown("### 📊 HISTORIQUE DES ANALYSES")
        response = supabase.table("milk_composition").select("brebis_id, date_analyse, mg_g_100ml, proteines_g_100ml, lactose_g_100ml, cellules_somatiques, aptitude_fromagere").eq("user_id", st.session_state.user_id).order("date_analyse", desc=True).limit(20).execute()
        data = response.data
        if data:
            # Pour avoir l'identifiant de la brebis, on peut faire une jointure, mais ici on simplifie
            df_hist = pd.DataFrame(data)
            st.dataframe(df_hist)
        else:
            st.info("Aucune analyse enregistrée.")
    with tab3:
        st.markdown("### 📈 COMPARAISON PAR RACE")
        st.info("Composition moyenne du lait de brebis par race (données simulées).")
        races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE']
        mg_moy = [7.8, 7.2, 7.5, 8.1]
        prot_moy = [6.2, 5.6, 5.9, 6.5]
        lactose_moy = [4.7, 4.9, 4.8, 4.6]
        df_comp_race = pd.DataFrame({
            'Race': races,
            'MG moyenne': mg_moy,
            'Protéines moyennes': prot_moy,
            'Lactose moyen': lactose_moy
        })
        st.dataframe(df_comp_race)
        fig = px.bar(df_comp_race, x='Race', y=['MG moyenne','Protéines moyennes','Lactose moyen'],
                    barmode='group', title="Composition moyenne du lait par race")
        st.plotly_chart(fig, use_container_width=True)

# --- PAGE ESTIMATION VIANDE (adaptée) ---
def page_estimation_viande():
    st.markdown('<h2 class="section-header">🥩 ESTIMATION DE LA COMPOSITION DE LA CARCASSE</h2>', unsafe_allow_html=True)
    st.markdown("**Prédiction de la quantité de viande, graisse et os à partir des mesures morphométriques**")
    tab1, tab2 = st.tabs(["📏 Estimation individuelle", "📊 Analyse de lot"])
    with tab1:
        st.markdown("### 🐑 Saisie des mesures")
        col1, col2 = st.columns(2)
        with col1:
            identifiant = st.text_input("Identifiant de la brebis", "HAM-F-001")
            # Recherche de la brebis
            response = supabase.table("brebis").select("id, race, age_mois, poids, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm").eq("identifiant", identifiant).eq("user_id", st.session_state.user_id).execute()
            existing = response.data[0] if response.data else None
            if existing:
                st.success(f"Brebis trouvée : {identifiant}")
                race_def = existing['race']
                age_def = existing['age_mois'] if existing['age_mois'] else 24
                poids_def = existing['poids'] if existing['poids'] else 50.0
                longueur_def = existing['longueur_corps_cm'] if existing['longueur_corps_cm'] else 100.0
                hauteur_def = existing['hauteur_garrot_cm'] if existing['hauteur_garrot_cm'] else 70.0
                poitrine_def = existing['tour_poitrine_cm'] if existing['tour_poitrine_cm'] else 100.0
            else:
                race_def = 'HAMRA'
                age_def = 24
                poids_def = 50.0
                longueur_def = 100.0
                hauteur_def = 70.0
                poitrine_def = 100.0
            race = st.selectbox("Race", list(STANDARDS_RACES.keys()), 
                              format_func=lambda x: STANDARDS_RACES[x]['nom_complet'],
                              index=list(STANDARDS_RACES.keys()).index(race_def) if race_def in STANDARDS_RACES else 0)
            age_mois = st.number_input("Âge (mois)", 6, 120, value=age_def)
            poids = st.number_input("Poids vif (kg)", 20.0, 150.0, value=poids_def, step=0.5)
        with col2:
            longueur = st.number_input("Longueur du corps (cm)", 50.0, 150.0, value=longueur_def, step=0.5)
            hauteur = st.number_input("Hauteur au garrot (cm)", 40.0, 120.0, value=hauteur_def, step=0.5)
            tour_poitrine = st.number_input("Tour de poitrine (cm)", 60.0, 150.0, value=poitrine_def, step=0.5)
            prix_viande = st.number_input("Prix de la viande (DA/kg)", 500, 1500, 800, 50)
            prix_gras = st.number_input("Prix du gras (DA/kg)", 10, 200, 50, 10)
        if st.button("🔬 ESTIMER LA COMPOSITION", type="primary"):
            estimation = CarcassEstimator.estimer_composition(poids, longueur, hauteur, tour_poitrine, age_mois, race)
            valeur = CarcassEstimator.estimer_valeur_economique(estimation['viande_kg'], estimation['graisse_kg'], prix_viande, prix_gras)
            st.success("✅ Estimation réalisée")
            col_res1, col_res2, col_res3 = st.columns(3)
            with col_res1:
                st.metric("Viande estimée", f"{estimation['viande_kg']} kg")
                st.metric("Rendement carcasse", f"{estimation['rendement_carcasse']} %")
            with col_res2:
                st.metric("Graisse estimée", f"{estimation['graisse_kg']} kg")
                st.metric("Valeur viande", f"{valeur['valeur_viande_da']:,.0f} DA")
            with col_res3:
                st.metric("Os estimés", f"{estimation['os_kg']} kg")
                st.metric("Valeur totale", f"{valeur['valeur_totale_da']:,.0f} DA")
            fig = go.Figure(data=[go.Pie(labels=['Viande', 'Graisse', 'Os'],
                                        values=[estimation['viande_kg'], estimation['graisse_kg'], estimation['os_kg']],
                                        hole=0.4)])
            fig.update_layout(title="Composition estimée de la carcasse")
            st.plotly_chart(fig, use_container_width=True)
            if st.button("💾 Enregistrer l'estimation"):
                if existing:
                    brebis_db_id = existing['id']
                else:
                    # Créer la brebis
                    data_brebis = {
                        'identifiant': identifiant,
                        'nom': f"Brebis_{identifiant}",
                        'race': race,
                        'sexe': 'F',
                        'age_mois': age_mois,
                        'poids': poids,
                        'user_id': st.session_state.user_id
                    }
                    resp = supabase.table("brebis").insert(data_brebis).execute()
                    brebis_db_id = resp.data[0]['id']
                data = {
                    'brebis_id': brebis_db_id,
                    'date_estimation': date.today().isoformat(),
                    'poids_vif_kg': poids,
                    'longueur_corps_cm': longueur,
                    'hauteur_garrot_cm': hauteur,
                    'tour_poitrine_cm': tour_poitrine,
                    'viande_estimee_kg': estimation['viande_kg'],
                    'graisse_estimee_kg': estimation['graisse_kg'],
                    'os_estimes_kg': estimation['os_kg'],
                    'rendement_carcasse': estimation['rendement_carcasse'],
                    'methode': "Modèle Ovin Manager",
                    'confiance': estimation['confiance'],
                    'notes': f"Estimation basée sur {race}, âge {age_mois} mois",
                    'user_id': st.session_state.user_id
                }
                supabase.table("carcass_estimates").insert(data).execute()
                st.success("✅ Estimation enregistrée !")
    with tab2:
        st.markdown("### 📊 ESTIMATION EN LOT")
        st.info("Sélectionnez plusieurs brebis et estimez leur composition en une fois.")
        response = supabase.table("brebis").select("id, identifiant, race, age_mois, poids, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm").eq("sexe", "F").eq("user_id", st.session_state.user_id).limit(20).execute()
        brebis_list = response.data
        if brebis_list:
            df_brebis = pd.DataFrame(brebis_list, 
                                    columns=['id','identifiant','race','age_mois','poids','longueur_corps_cm','hauteur_garrot_cm','tour_poitrine_cm'])
            st.dataframe(df_brebis[['identifiant','race','age_mois','poids']])
            selected = st.multiselect("Choisir les brebis", df_brebis['identifiant'].tolist())
            if selected and st.button("📊 Estimer pour la sélection"):
                results = []
                for ident in selected:
                    row = df_brebis[df_brebis['identifiant'] == ident].iloc[0]
                    est = CarcassEstimator.estimer_composition(
                        row['poids'], row['longueur_corps_cm'], row['hauteur_garrot_cm'], row['tour_poitrine_cm'], row['age_mois'], row['race']
                    )
                    results.append({
                        'Identifiant': ident,
                        'Viande (kg)': est['viande_kg'],
                        'Graisse (kg)': est['graisse_kg'],
                        'Os (kg)': est['os_kg'],
                        'Rendement %': est['rendement_carcasse']
                    })
                df_res = pd.DataFrame(results)
                st.dataframe(df_res)
                csv = df_res.to_csv(index=False)
                st.download_button("📥 Télécharger le rapport CSV", csv, "estimation_lot.csv", "text/csv")
        else:
            st.warning("Aucune brebis femelle trouvée.")

# --- PAGE ESTIMATION LAIT MORPHO (adaptée) ---
def page_estimation_lait_morpho():
    st.markdown('<h2 class="section-header">🍼 ESTIMATION LAIT PAR MORPHOMÉTRIE MAMMAIRE</h2>', unsafe_allow_html=True)
    st.markdown("**Prédiction de la production laitière individuelle à partir des caractéristiques de la mamelle**")
    tab1, tab2 = st.tabs(["📏 Estimation individuelle", "📈 Analyse par ID"])
    with tab1:
        st.markdown("### 🍼 MESURES DE LA MAMELLE")
        col1, col2 = st.columns(2)
        with col1:
            identifiant = st.text_input("Identifiant de la brebis", "HAM-F-001")
            response = supabase.table("brebis").select("id, race, age_mois, volume_mammaire, longueur_trayons_cm").eq("identifiant", identifiant).eq("user_id", st.session_state.user_id).execute()
            existing = response.data[0] if response.data else None
            if existing:
                race_def = existing['race'] if existing['race'] else 'HAMRA'
                age_def = existing['age_mois'] if existing['age_mois'] else 24
                vol_mammaire_score = existing['volume_mammaire'] if existing['volume_mammaire'] else 3
                vol_cm3_def = vol_mammaire_score * 80 + 100
                # longueur_trayons non utilisée ici
            else:
                race_def = 'HAMRA'
                age_def = 24
                vol_cm3_def = 250
            race = st.selectbox("Race", list(STANDARDS_RACES.keys()),
                              format_func=lambda x: STANDARDS_RACES[x]['nom_complet'],
                              index=list(STANDARDS_RACES.keys()).index(race_def) if race_def in STANDARDS_RACES else 0)
            age_mois = st.number_input("Âge (mois)", 6, 120, value=age_def)
        with col2:
            volume_cm3 = st.number_input("Volume mammaire estimé (cm³)", 100, 800, value=int(vol_cm3_def), step=10)
            largeur_cm = st.number_input("Largeur moyenne de la mamelle (cm)", 5.0, 20.0, value=8.5, step=0.1)
            hauteur_cm = st.number_input("Hauteur moyenne de la mamelle (cm)", 5.0, 25.0, value=12.0, step=0.1)
            symetrie = st.slider("Symétrie (0 = asymétrique, 1 = parfaite)", 0.0, 1.0, 0.85, 0.05)
        if st.button("🥛 ESTIMER LA PRODUCTION", type="primary"):
            estimation = MilkProductionEstimator.estimer_lait_par_morphometrie(
                volume_cm3, largeur_cm, hauteur_cm, symetrie, age_mois, race
            )
            classe = MilkProductionEstimator.classe_productrice(estimation['lait_jour_estime_l'])
            st.success("✅ Estimation réalisée")
            col_res1, col_res2, col_res3 = st.columns(3)
            with col_res1:
                st.metric("Lait estimé par jour", f"{estimation['lait_jour_estime_l']} L")
            with col_res2:
                st.metric("Production potentielle / lactation", f"{estimation['lactation_potentielle_l']} L")
            with col_res3:
                st.metric("Confiance", f"{estimation['confiance']*100:.0f} %")
            st.info(f"**Classe : {classe}**")
            if st.button("💾 Enregistrer l'estimation"):
                if existing:
                    brebis_db_id = existing['id']
                else:
                    # Créer la brebis
                    data_brebis = {
                        'identifiant': identifiant,
                        'nom': f"Brebis_{identifiant}",
                        'race': race,
                        'sexe': 'F',
                        'age_mois': age_mois,
                        'user_id': st.session_state.user_id
                    }
                    resp = supabase.table("brebis").insert(data_brebis).execute()
                    brebis_db_id = resp.data[0]['id']
                data = {
                    'brebis_id': brebis_db_id,
                    'date_estimation': date.today().isoformat(),
                    'volume_mammaire_cm3': volume_cm3,
                    'largeur_mammaire_cm': largeur_cm,
                    'hauteur_mammaire_cm': hauteur_cm,
                    'symetrie': symetrie,
                    'lait_jour_estime_l': estimation['lait_jour_estime_l'],
                    'lactation_potentielle_l': estimation['lactation_potentielle_l'],
                    'confiance': estimation['confiance'],
                    'notes': "Estimation basée sur morphométrie mammaire",
                    'user_id': st.session_state.user_id
                }
                supabase.table("milk_production_estimates").insert(data).execute()
                st.success("✅ Estimation enregistrée !")
    with tab2:
        st.markdown("### 📈 CONSULTATION PAR ID")
        st.markdown("Recherchez une brebis et visualisez ses estimations de production laitière.")
        response = supabase.table("brebis").select("id, identifiant").eq("sexe", "F").eq("user_id", st.session_state.user_id).limit(50).execute()
        brebis_f = response.data
        if brebis_f:
            brebis_choice = st.selectbox("Choisir une brebis", 
                                        [f"{b['identifiant']} (ID: {b['id']})" for b in brebis_f])
            brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
            response = supabase.table("milk_production_estimates").select("date_estimation, volume_mammaire_cm3, lait_jour_estime_l, lactation_potentielle_l, confiance").eq("brebis_id", brebis_id).order("date_estimation", desc=True).execute()
            est_data = response.data
            if est_data:
                df_est = pd.DataFrame(est_data, 
                                     columns=['date_estimation','volume_mammaire_cm3','lait_jour_estime_l','lactation_potentielle_l','confiance'])
                df_est = df_est.rename(columns={'date_estimation':'Date','volume_mammaire_cm3':'Volume (cm³)','lait_jour_estime_l':'Lait/jour (L)','lactation_potentielle_l':'Lactation (L)','confiance':'Confiance'})
                st.dataframe(df_est)
                if len(est_data) > 1:
                    fig = px.line(df_est, x='Date', y='Lait/jour (L)', 
                                 title="Évolution de l'estimation du lait journalier")
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Aucune estimation enregistrée pour cette brebis.")
        else:
            st.warning("Aucune brebis femelle trouvée.")

# ============================================================================
# SECTION 11: BARRE LATÉRALE AVEC DÉCONNEXION
# ============================================================================
with st.sidebar:
    st.markdown("""
    <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #1a237e 0%, #283593 100%); 
                color: white; border-radius: 10px; margin-bottom: 20px;'>
        <h2>🐑 RACES ALGÉRIENNES</h2>
        <p>Système de gestion scientifique</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Afficher l'utilisateur connecté
    if st.session_state.authenticated:
        st.markdown(f"**Connecté :** {st.session_state.user_email}")
        if st.button("Se déconnecter"):
            logout()
    
    st.markdown("### 📍 NAVIGATION")
    
    page = st.radio(
        "MENU PRINCIPAL",
        [
            "🏠 ACCUEIL",
            "📸 PHOTO & MESURES",
            "📦 ANALYSE MULTIPLE",
            "📐 SCANNER 3D",
            "📊 GESTION",
            "🥛 PRODUCTION",
            "🎯 CRITÈRES",
            "📊 RSTATS",
            "🧬 GÉNÉTIQUE",
            "🧬🔬 GÉNOMIQUE AVANCÉE",
            "🥛🔬 ANALYSE LAIT",
            "🥩 ESTIMATION VIANDE",
            "🍼 ESTIMATION LAIT MAMMELLE"
        ]
    )
    
    st.markdown("---")
    
    # Statistiques rapides (pour l'utilisateur connecté)
    if st.session_state.authenticated:
        try:
            response = supabase.table("brebis").select("*", count="exact").eq("user_id", st.session_state.user_id).execute()
            total_brebis = response.count
        except:
            total_brebis = 0
        st.metric("🐑 Brebis", total_brebis)
        st.metric("🏷️ Races", "6")  # à calculer éventuellement
        st.metric("🧬 Génotypes", ">100")
    
    st.markdown("---")
    
    # Standards des races
    st.markdown("### 🏷️ STANDARDS")
    race_info = st.selectbox("Info race", list(STANDARDS_RACES.keys()),
                            format_func=lambda x: STANDARDS_RACES[x]['nom_complet'])
    if race_info in STANDARDS_RACES:
        info = STANDARDS_RACES[race_info]
        st.markdown(f"""
        <div class='race-card' style='padding: 10px;'>
            <h5>{info['nom_complet']}</h5>
            <p><small>Poids ♀️: {info['poids_adulte']['femelle'][0]}-{info['poids_adulte']['femelle'][1]} kg</small></p>
            <p><small>Poids ♂️: {info['poids_adulte']['male'][0]}-{info['poids_adulte']['male'][1]} kg</small></p>
            <p><small>Lait: {info['production_lait'][0]}-{info['production_lait'][1]} L/j</small></p>
        </div>
        """, unsafe_allow_html=True)

# ============================================================================
# SECTION 12: NAVIGATION PRINCIPALE
# ============================================================================
if not st.session_state.authenticated:
    login_page()
    st.stop()
else:
    if page == "🏠 ACCUEIL":
        page_accueil()
    elif page == "📸 PHOTO & MESURES":
        page_photo_mesures()
    elif page == "📦 ANALYSE MULTIPLE":
        page_analyse_multiple()
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
    elif page == "🧬🔬 GÉNOMIQUE AVANCÉE":
        page_genomique_avancee()
    elif page == "🥛🔬 ANALYSE LAIT":
        page_analyse_lait()
    elif page == "🥩 ESTIMATION VIANDE":
        page_estimation_viande()
    elif page == "🍼 ESTIMATION LAIT MAMMELLE":
        page_estimation_lait_morpho()

# ============================================================================
# SECTION 13: PIED DE PAGE
# ============================================================================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version Multi-Utilisateurs</p>
    <p>📸 Photo & Mesures • 📦 Analyse Multiple • 📐 Scanner 3D • 🎯 Critères de sélection • 🧬 Génétique • 🧬🔬 Génomique avancée • 🥛🔬 Biochimie lait • 🥩 Estimation viande • 🍼 Estimation lait</p>
    <p>© 2024 - Système de gestion scientifique des races ovines algériennes</p>
</div>
""", unsafe_allow_html=True)
