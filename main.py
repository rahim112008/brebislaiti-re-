"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D et Génétique
Base de données simulée de races ovines algériennes
Version avec critères de sélection mammaires et noms génériques
CODE COMPLET AVEC MODULE PHOTO & MESURES AUTOMATIQUES
"""

# ============================================================================
# SECTION 1: IMPORTS - AJOUTER IMPORT cv2
# ============================================================================
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
import random
import math
import io
import base64
from PIL import Image, ImageDraw
import tempfile
import os
import cv2 

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
    .photo-card {
        background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(30,60,114,0.2);
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# SECTION 4: STANDARDS DES RACES ALGÉRIENNES - VERSION SÉCURISÉE
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

# Fonction helper pour accéder aux données de manière sécurisée
def get_race_data(race, key, default=None):
    """Récupère les données d'une race de manière sécurisée"""
    if race in STANDARDS_RACES:
        data = STANDARDS_RACES[race]
        if key in data:
            return data[key]
    
    # Retourne des valeurs par défaut si la race n'existe pas
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
# SECTION 5: MODULE PHOTO & MESURES - VERSION AVEC 2 PHOTOS
# ============================================================================

class OvinPhotoAnalyzer:
    """Analyseur pour photos profil ET arrière"""
    
    def __init__(self):
        self.etalon_type = None
        self.etalon_size_mm = None
        self.pixel_per_mm = None
        self.profile_image = None
        self.rear_image = None
        
    def set_etalon(self, etalon_type):
        """Définit l'étalon de référence"""
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
        """Détecte l'étalon dans n'importe quelle photo"""
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
                box = np.int0(box)
                
                width = np.linalg.norm(box[0] - box[1])
                height = np.linalg.norm(box[1] - box[2])
                longest_side = max(width, height)
                
                self.pixel_per_mm = longest_side / self.etalon_size_mm
                return {'type': 'rectangle', 'box': box}
            
            return None
        except Exception as e:
            st.warning(f"Détection étalon: {str(e)}")
            return None
    
    # === MESURES DU CORPS (PHOTO DE PROFIL) ===
    def analyze_profile_photo(self, image):
        """Analyse la photo de profil pour mesures corporelles"""
        measurements = {}
        
        if self.pixel_per_mm is None:
            etalon_info = self.detect_etalon(image)
            if not etalon_info:
                return None
        
        try:
            # Détection du corps
            gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
            thresh = cv2.adaptiveThreshold(gray, 255, 
                                          cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                          cv2.THRESH_BINARY_INV, 11, 2)
            
            contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            if contours:
                animal_contour = max(contours, key=cv2.contourArea)
                x, y, w, h = cv2.boundingRect(animal_contour)
                
                # 1. Longueur du corps
                measurements['longueur_corps_cm'] = (w / self.pixel_per_mm) / 10
                
                # 2. Hauteur au garrot
                measurements['hauteur_garrot_cm'] = (h / self.pixel_per_mm) / 10
                
                # 3. Estimation tour de poitrine
                middle_y = y + h // 2
                chest_width = 0
                for i in range(x, x + w):
                    column_pixels = []
                    for j in range(y, y + h):
                        if thresh[j, i] > 0:
                            column_pixels.append(j)
                    if column_pixels:
                        col_height = max(column_pixels) - min(column_pixels)
                        if col_height > chest_width:
                            chest_width = col_height
                
                measurements['tour_poitrine_cm'] = (chest_width * np.pi / self.pixel_per_mm) / 10
                
                # 4. Indices corporels
                measurements['ratio_longueur_hauteur'] = measurements['longueur_corps_cm'] / measurements['hauteur_garrot_cm']
                
                return measurements
        
        except Exception as e:
            st.warning(f"Analyse profil échouée: {str(e)}")
        
        return None
    
    # === MESURES DES MAMELLES (PHOTO ARRIÈRE) ===
    def analyze_rear_photo(self, image, is_female=True):
        """Analyse la photo arrière pour évaluer les mamelles"""
        if not is_female:
            return {'message': 'Animal mâle - pas de mamelles à évaluer'}
        
        mammary_data = {}
        
        if self.pixel_per_mm is None:
            etalon_info = self.detect_etalon(image)
            if not etalon_info:
                return {'error': 'Étalon non détecté'}
        
        try:
            # Prétraitement
            gray = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
            blurred = cv2.GaussianBlur(gray, (5, 5), 0)
            
            # Seuillage pour détecter les mamelles (plus claires que le corps)
            _, thresh = cv2.threshold(blurred, 150, 255, cv2.THRESH_BINARY)
            
            # Recherche dans le tiers inférieur (région mammaire)
            height, width = gray.shape
            search_region = thresh[2*height//3:height, width//4:3*width//4]
            
            # Détection des contours
            contours, _ = cv2.findContours(search_region, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            mammary_contours = []
            for contour in contours:
                area = cv2.contourArea(contour)
                # Filtre de taille pour les mamelles
                if 200 < area < 3000:
                    mammary_contours.append(contour)
            
            mammary_data['nombre_mamelles_detectees'] = len(mammary_contours)
            
            if mammary_contours:
                # Mesurer chaque mamelle
                volumes = []
                widths = []
                heights = []
                
                for i, contour in enumerate(mammary_contours[:4]):  # Max 4 mamelles
                    if len(contour) >= 5:
                        ellipse = cv2.fitEllipse(contour)
                        (x_ell, y_ell), (w_ell, h_ell), angle = ellipse
                        
                        # Convertir en mm
                        width_mm = w_ell / self.pixel_per_mm
                        height_mm = h_ell / self.pixel_per_mm
                        
                        # Volume estimé (ellipsoïde)
                        volume = (4/3) * np.pi * (width_mm/2) * (height_mm/2) * (height_mm/2)
                        
                        mammary_data[f'mamelle_{i+1}_largeur_mm'] = width_mm
                        mammary_data[f'mamelle_{i+1}_hauteur_mm'] = height_mm
                        mammary_data[f'mamelle_{i+1}_volume_mm3'] = volume
                        
                        volumes.append(volume)
                        widths.append(width_mm)
                        heights.append(height_mm)
                
                if volumes:
                    # Calcul des moyennes
                    mammary_data['volume_mammaire_moyen_cm3'] = np.mean(volumes) / 1000
                    mammary_data['largeur_mammaire_moyenne_cm'] = np.mean(widths) / 10
                    mammary_data['hauteur_mammaire_moyenne_cm'] = np.mean(heights) / 10
                    
                    # Calcul de la symétrie (si 2 mamelles ou plus)
                    if len(volumes) >= 2:
                        left_volumes = volumes[:len(volumes)//2]
                        right_volumes = volumes[len(volumes)//2:]
                        
                        if left_volumes and right_volumes:
                            avg_left = np.mean(left_volumes)
                            avg_right = np.mean(right_volumes)
                            symmetry = min(avg_left, avg_right) / max(avg_left, avg_right) if max(avg_left, avg_right) > 0 else 0
                            mammary_data['symetrie_mammaire'] = round(symmetry, 2)
                    
                    # Score de développement mammaire
                    mammary_data['score_developpement'] = self.calculate_mammary_score(volumes, widths, heights)
            
            return mammary_data
            
        except Exception as e:
            return {'error': f'Analyse mamelles échouée: {str(e)}'}
    
    def calculate_mammary_score(self, volumes, widths, heights):
        """Calcule un score de développement mammaire (1-10)"""
        if not volumes:
            return 0
        
        avg_volume = np.mean(volumes) / 1000  # en cm³
        avg_width = np.mean(widths) / 10      # en cm
        avg_height = np.mean(heights) / 10    # en cm
        
        # Critères de scoring (basés sur standards ovins)
        volume_score = min(10, avg_volume / 50)  # 500 cm³ = score 10
        width_score = min(10, avg_width / 1.5)   # 15 cm = score 10
        height_score = min(10, avg_height / 2.0) # 20 cm = score 10
        
        return round((volume_score + width_score + height_score) / 3, 1)
    
    def get_mammary_classification(self, mammary_data):
        """Classifie l'aptitude laitière basée sur les mamelles"""
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

# ============================================================================
# MODIFICATION DE LA PAGE PHOTO_MESURES
# ============================================================================

def page_photo_mesures():
    """Page de capture photo avec 2 vues (profil + arrière)"""
    st.markdown('<h2 class="section-header">📸 CARACTÉRISATION COMPLÈTE DES BREBIS LAITIÈRES</h2>', unsafe_allow_html=True)
    
    # Initialisation
    if 'photo_analyzer' not in st.session_state:
        st.session_state.photo_analyzer = OvinPhotoAnalyzer()
    
    # Onglets pour les 2 photos
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
                }[x]
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
            </div>
            """, unsafe_allow_html=True)
        
        # Mettre à jour l'analyseur
        st.session_state.photo_analyzer.set_etalon(etalon_type)
        
        # Instructions visuelles
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
            
            # Schéma ASCII simplifié
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
        
        # Capture photo de profil
        profile_img = st.camera_input(
            "Prenez une photo de PROFIL de la brebis",
            key="camera_profile",
            help="Photo latérale pour mesures du corps"
        )
        
        if profile_img is not None:
            # Sauvegarder l'image
            bytes_data = profile_img.getvalue()
            image = Image.open(io.BytesIO(bytes_data))
            st.session_state.profile_image = np.array(image)
            
            st.success("✅ Photo de profil enregistrée!")
            st.image(image, caption="Photo de profil - Mesures corporelles", use_column_width=True)
            
            # Bouton d'analyse immédiate
            if st.button("📏 Analyser les mesures corporelles", type="primary"):
                with st.spinner("Analyse en cours..."):
                    measurements = st.session_state.photo_analyzer.analyze_profile_photo(
                        st.session_state.profile_image
                    )
                    
                    if measurements:
                        st.session_state.body_measurements = measurements
                        st.success("✅ Mesures corporelles extraites!")
                        
                        # Afficher les résultats
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
                        
                        # Enregistrer dans session
                        st.session_state.has_profile_analysis = True
                    else:
                        st.error("❌ Impossible d'analyser la photo. Vérifiez l'étalon.")
    
    with tab3:
        st.markdown("### 🍼 2ÈME PHOTO : VUE ARRIÈRE")
        st.markdown("*Pour l'évaluation des mamelles (volume, symétrie, développement)*")
        
        # Vérifier si on a déjà la photo de profil
        if 'has_profile_analysis' not in st.session_state:
            st.warning("⚠️ Prenez d'abord la photo de profil (onglet 2)")
            return
        
        # Capture photo arrière
        rear_img = st.camera_input(
            "Prenez une photo ARRIÈRE de la brebis",
            key="camera_rear",
            help="Photo postérieure pour évaluer les mamelles"
        )
        
        if rear_img is not None:
            # Sauvegarder l'image
            bytes_data = rear_img.getvalue()
            image = Image.open(io.BytesIO(bytes_data))
            st.session_state.rear_image = np.array(image)
            
            st.success("✅ Photo arrière enregistrée!")
            st.image(image, caption="Photo arrière - Évaluation des mamelles", use_column_width=True)
            
            # Formulaire pour informations
            with st.form("mammary_info_form"):
                st.markdown("#### ℹ️ INFORMATIONS COMPLÉMENTAIRES")
                
                col_info1, col_info2 = st.columns(2)
                
                with col_info1:
                    race = st.selectbox("Race de la brebis", 
                                       list(STANDARDS_RACES.keys()),
                                       format_func=lambda x: STANDARDS_RACES[x]['nom_complet'])
                    
                    age_mois = st.number_input("Âge (mois)", 6, 120, 24)
                
                with col_info2:
                    poids = st.number_input("Poids estimé (kg)", 20.0, 100.0, 50.0, 0.5)
                    
                    # Demander si c'est une femelle
                    is_female = st.radio("Sexe", ["Femelle", "Mâle"], horizontal=True)
                
                analyse_mammary = st.form_submit_button("🔍 ANALYSER LES MAMELLES", type="primary")
                
                if analyse_mammary:
                    with st.spinner("Analyse détaillée des mamelles..."):
                        # Analyser les mamelles
                        mammary_data = st.session_state.photo_analyzer.analyze_rear_photo(
                            st.session_state.rear_image,
                            is_female=(is_female == "Femelle")
                        )
                        
                        if 'error' not in mammary_data:
                            st.session_state.mammary_data = mammary_data
                            
                            # AFFICHAGE DES RÉSULTATS COMPLETS
                            st.markdown("---")
                            st.markdown("## 📊 RAPPORT COMPLET DE CARACTÉRISATION")
                            
                            # Section 1: Corps
                            st.markdown("### 🐑 CARACTÉRISTIQUES CORPORELES")
                            if 'body_measurements' in st.session_state:
                                body_df = pd.DataFrame([
                                    {"Paramètre": "Longueur corps", "Valeur": f"{st.session_state.body_measurements['longueur_corps_cm']:.1f} cm", "Norme": "90-130 cm"},
                                    {"Paramètre": "Hauteur garrot", "Valeur": f"{st.session_state.body_measurements['hauteur_garrot_cm']:.1f} cm", "Norme": "60-90 cm"},
                                    {"Paramètre": "Tour poitrine", "Valeur": f"{st.session_state.body_measurements['tour_poitrine_cm']:.1f} cm", "Norme": "95-130 cm"},
                                    {"Paramètre": "Ratio L/H", "Valeur": f"{st.session_state.body_measurements['ratio_longueur_hauteur']:.2f}", "Norme": "1.4-1.6"}
                                ])
                                st.dataframe(body_df, use_container_width=True, hide_index=True)
                            
                            # Section 2: Mamelles
                            st.markdown("### 🍼 CARACTÉRISTIQUES MAMMAIRES")
                            
                            if mammary_data.get('nombre_mamelles_detectees', 0) > 0:
                                # Afficher les mesures détaillées
                                mammary_metrics = []
                                
                                if 'volume_mammaire_moyen_cm3' in mammary_data:
                                    mammary_metrics.append({
                                        "Paramètre": "Volume moyen",
                                        "Valeur": f"{mammary_data['volume_mammaire_moyen_cm3']:.1f} cm³",
                                        "Évaluation": self.evaluate_volume(mammary_data['volume_mammaire_moyen_cm3'])
                                    })
                                
                                if 'largeur_mammaire_moyenne_cm' in mammary_data:
                                    mammary_metrics.append({
                                        "Paramètre": "Largeur moyenne",
                                        "Valeur": f"{mammary_data['largeur_mammaire_moyenne_cm']:.1f} cm",
                                        "Évaluation": self.evaluate_width(mammary_data['largeur_mammaire_moyenne_cm'])
                                    })
                                
                                if 'symetrie_mammaire' in mammary_data:
                                    mammary_metrics.append({
                                        "Paramètre": "Symétrie",
                                        "Valeur": f"{mammary_data['symetrie_mammaire']:.2f}",
                                        "Évaluation": self.evaluate_symmetry(mammary_data['symetrie_mammaire'])
                                    })
                                
                                if mammary_metrics:
                                    mammary_df = pd.DataFrame(mammary_metrics)
                                    st.dataframe(mammary_df, use_container_width=True, hide_index=True)
                                
                                # Classification finale
                                classification = st.session_state.photo_analyzer.get_mammary_classification(mammary_data)
                                
                                st.markdown("#### 🎯 APTITUDE LAITIÈRE")
                                
                                # Afficher avec couleur selon le score
                                if "EXCELLENT" in classification:
                                    st.success(f"**{classification}**")
                                elif "BON" in classification:
                                    st.info(f"**{classification}**")
                                elif "MOYEN" in classification:
                                    st.warning(f"**{classification}**")
                                else:
                                    st.error(f"**{classification}**")
                                
                                # Recommandations
                                st.markdown("#### 💡 RECOMMANDATIONS")
                                
                                if mammary_data.get('score_developpement', 0) >= 6:
                                    st.info("""
                                    ✅ **Bonne candidate pour :**
                                    - Reproduction
                                    - Production laitière
                                    - Amélioration génétique
                                    """)
                                else:
                                    st.warning("""
                                    ⚠️ **À surveiller ou réformer :**
                                    - Contrôler régulièrement
                                    - Évaluer la production réelle
                                    - Envisager le renouvellement
                                    """)
                            
                            else:
                                st.warning("Aucune mamelle détectée. Vérifiez la position de l'animal.")
                            
                            # Bouton d'enregistrement
                            st.markdown("---")
                            if st.button("💾 ENREGISTRER DANS LA BASE DE DONNÉES", type="primary", use_container_width=True):
                                self.save_complete_characterization(race, age_mois, poids, is_female)
                        
                        else:
                            st.error(f"Erreur: {mammary_data['error']}")
    
    # Fonctions d'évaluation
    def evaluate_volume(self, volume_cm3):
        if volume_cm3 > 400: return "Très bon"
        elif volume_cm3 > 250: return "Bon"
        elif volume_cm3 > 150: return "Moyen"
        else: return "Faible"
    
    def evaluate_width(self, width_cm):
        if width_cm > 12: return "Large"
        elif width_cm > 8: return "Normale"
        else: return "Étroite"
    
    def evaluate_symmetry(self, symmetry):
        if symmetry > 0.9: return "Excellente"
        elif symmetry > 0.8: return "Bonne"
        elif symmetry > 0.7: return "Acceptable"
        else: return "Asymétrique"
    
    def save_complete_characterization(self, race, age_mois, poids, sexe):
        """Sauvegarde toutes les données dans la base"""
        try:
            cursor = conn.cursor()
            
            # Générer un ID unique
            identifiant = f"{race[:3]}-{sexe[:1]}-{datetime.now().strftime('%y%m%d')}"
            
            # Préparer les données
            complete_data = {
                'identifiant': identifiant,
                'race': race,
                'age_mois': age_mois,
                'poids': poids,
                'sexe': 'F' if sexe == "Femelle" else 'M',
                'body_measurements': json.dumps(st.session_state.body_measurements) if 'body_measurements' in st.session_state else None,
                'mammary_data': json.dumps(st.session_state.mammary_data) if 'mammary_data' in st.session_state else None,
                'date_analyse': datetime.now().isoformat()
            }
            
            # Insérer dans la base
            cursor.execute('''
                INSERT INTO brebis (
                    identifiant, nom, race, sexe, age_mois, poids,
                    longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm,
                    volume_mammaire, symetrie_mammaire, notes, statut
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                identifiant,
                f"Char_{identifiant}",
                race,
                'F' if sexe == "Femelle" else 'M',
                age_mois,
                poids,
                st.session_state.body_measurements.get('longueur_corps_cm', 0) if 'body_measurements' in st.session_state else 0,
                st.session_state.body_measurements.get('hauteur_garrot_cm', 0) if 'body_measurements' in st.session_state else 0,
                st.session_state.body_measurements.get('tour_poitrine_cm', 0) if 'body_measurements' in st.session_state else 0,
                st.session_state.mammary_data.get('score_developpement', 0) if 'mammary_data' in st.session_state else 0,
                st.session_state.mammary_data.get('symetrie_mammaire', 0) if 'mammary_data' in st.session_state else 0,
                f"Caractérisation complète - {classification if 'classification' in locals() else 'Non classé'}",
                'active'
            ))
            
            conn.commit()
            st.success(f"✅ Données enregistrées pour {identifiant}!")
            
            # Option de téléchargement
            st.download_button(
                label="📥 Télécharger le rapport complet (JSON)",
                data=json.dumps(complete_data, indent=2, ensure_ascii=False),
                file_name=f"caracterisation_{identifiant}.json",
                mime="application/json"
            )
            
        except Exception as e:
            st.error(f"❌ Erreur d'enregistrement: {str(e)}")

# ============================================================================
# SECTION 6: FONCTIONS STATISTIQUES (sans scipy)
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
# SECTION 7: BASE DE DONNÉES - VERSION SÉCURISÉE
# ============================================================================
def init_database_safe():
    """Initialise la base de données avec gestion robuste des erreurs"""
    try:
        # Utiliser un fichier temporaire pour Streamlit Cloud
        temp_db = tempfile.NamedTemporaryFile(delete=False, suffix='.db')
        db_path = temp_db.name
        temp_db.close()
        
        conn = sqlite3.connect(db_path, check_same_thread=False)
        cursor = conn.cursor()
        
        # Table brebis avec critères mammaires
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS brebis (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                identifiant TEXT UNIQUE NOT NULL,
                nom TEXT,
                race TEXT,
                sous_race TEXT,
                sexe TEXT,
                date_naissance DATE,
                age_mois INTEGER,
                poids FLOAT,
                score_condition INTEGER,
                
                -- Caractères morphologiques
                couleur_robe TEXT,
                intensite_couleur INTEGER,
                cornes BOOLEAN,
                taille_cornes_cm FLOAT,
                forme_cornes TEXT,
                type_laine TEXT,
                qualite_laine INTEGER,
                
                -- Mensurations
                longueur_corps_cm FLOAT,
                hauteur_garrot_cm FLOAT,
                largeur_bassin_cm FLOAT,
                tour_poitrine_cm FLOAT,
                circonference_tete_cm FLOAT,
                longueur_oreille_cm FLOAT,
                
                -- Critères mammaires (nouveau)
                volume_mammaire INTEGER,
                symetrie_mammaire INTEGER,
                insertion_trayons INTEGER,
                longueur_trayons_cm FLOAT,
                orientation_trayons TEXT,
                
                -- Caractères qualitatifs
                temperement TEXT,
                aptitude TEXT,
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
                notes TEXT
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
                notes TEXT
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
                date_analyse DATE
            )
        ''')
        
        # Vérifier si la base est vide et la peupler
        cursor.execute("SELECT COUNT(*) FROM brebis")
        count = cursor.fetchone()[0]
        
        if count == 0:
            peupler_base_races_safe(cursor, conn)
        
        conn.commit()
        return conn
        
    except Exception as e:
        st.error(f"Erreur d'initialisation: {str(e)}")
        # Retourner une connexion mémoire en cas d'erreur
        conn = sqlite3.connect(':memory:', check_same_thread=False)
        return conn

def peupler_base_races_safe(cursor, conn):
    """Peuple la base avec des races algériennes - VERSION SÉCURISÉE"""
    races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU']
    brebis_data = []
    
    for i in range(1, 21):  # Réduit à 20 animaux pour plus de rapidité
        race = random.choice(races)
        sexe = random.choice(['F', 'M'])
        
        # Générer identifiant et nom générique
        race_code = race[:3] if race != 'INCONNU' else 'INC'
        identifiant = f"{race_code}-{sexe}-2023-{i:03d}"
        
        # Noms génériques
        if sexe == 'F':
            nom = f"F{race_code}{i:03d}"
        else:
            nom = f"M{race_code}{i:03d}"
        
        age_mois = random.randint(12, 84)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        
        # Poids selon race et sexe - VERSION SÉCURISÉE
        try:
            poids_data = get_race_data(race, 'poids_adulte')
            if sexe.lower() in poids_data:
                poids_min, poids_max = poids_data[sexe.lower()]
            else:
                poids_min, poids_max = (35, 60) if sexe == 'F' else (50, 80)
        except:
            poids_min, poids_max = (35, 60) if sexe == 'F' else (50, 80)
        
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
        couleur_robe = random.choice(couleurs.get(race, ['Indéterminée']))
        
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
        
        # Score conformation
        score_conformation = random.uniform(5.0, 9.0)
        
        # Mensurations - VERSION SÉCURISÉE
        try:
            mensurations = get_race_data(race, 'mensurations')
            longueur_corps = random.uniform(*mensurations['longueur_cm'])
            hauteur_garrot = random.uniform(*mensurations['hauteur_cm'])
            largeur_bassin = random.uniform(*mensurations['largeur_bassin_cm'])
            tour_poitrine = random.uniform(*mensurations['tour_poitrine_cm'])
        except:
            longueur_corps = random.uniform(80, 120)
            hauteur_garrot = random.uniform(55, 80)
            largeur_bassin = random.uniform(30, 50)
            tour_poitrine = random.uniform(85, 120)
        
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
    try:
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
        
    except Exception as e:
        # Insérer une seule donnée de test
        cursor.execute('''
            INSERT INTO brebis (identifiant, nom, race, sexe, age_mois, poids, statut)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', ('TEST-001', 'TestBrebis', 'HAMRA', 'F', 24, 55.5, 'active'))
        conn.commit()

# Initialiser la base de données
@st.cache_resource
def get_database_connection():
    """Obtient une connexion à la base de données avec cache"""
    return init_database_safe()

# Obtenir la connexion
conn = get_database_connection()

# ============================================================================
# SECTION 8: MODULE SCANNER 3D
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
            'HAMRA': (139, 0, 0),
            'OUDA': (255, 255, 255),
            'SIDAHOU': (50, 50, 50),
            'BERBERE': (165, 42, 42),
            'CROISE': (160, 120, 80),
            'INCONNU': (200, 200, 200)
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
        
        # Informations
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {race}", fill='black')
        draw.text((10, 50), f"Poids: {brebis_info.get('poids', 0):.1f} kg", fill='black')
        
        return image
    
    @staticmethod
    def simuler_scan_3d(brebis_info):
        """Simule un scan 3D réaliste"""
        np.random.seed(hash(str(brebis_info.get('identifiant', ''))) % 10000)
        
        n_points = 200  # Réduit pour la performance
        points = []
        
        poids = brebis_info.get('poids', 50)
        
        # Rayons approximatifs
        rx = 0.6 * poids**0.33
        ry = 1.2 * poids**0.33
        rz = 0.8 * poids**0.33
        
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
# SECTION 9: MODULE GÉNÉTIQUE
# ============================================================================
class ModuleGenetique:
    """Module d'analyse génétique"""
    
    @staticmethod
    def generer_genotype(brebis_id, race):
        """Génère un génotype simulé"""
        genotypes = []
        
        for i in range(5):  # Réduit pour la performance
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
        """Calcule la diversité génétique"""
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
# SECTION 10: PAGE ACCUEIL
# ============================================================================
def page_accueil():
    """Page d'accueil avec vue d'ensemble"""
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - RACES ALGÉRIENNES</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique des races ovines algériennes**")
    
    # Métriques principales
    try:
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
    
    except Exception as e:
        st.info("Bienvenue dans Ovin Manager Pro! Le système est en cours d'initialisation.")
        st.markdown("### Fonctionnalités disponibles:")
        st.markdown("""
        - **📸 Photo & Mesures**: Capture photo avec étalon et mesures automatiques
        - **📐 Scanner 3D**: Simulation de scans 3D
        - **📊 Gestion**: Suivi du troupeau
        - **🥛 Production**: Suivi laitier
        - **🎯 Critères**: Évaluation des mamelles
        - **📊 Statistiques**: Analyses avancées
        - **🧬 Génétique**: Analyses génomiques
        """)

# ============================================================================
# SECTION 11: PAGE SCANNER 3D
# ============================================================================
def page_scanner_3d():
    """Page du scanner 3D avec saisie manuelle"""
    st.markdown('<h2 class="section-header">📐 SCANNER 3D & SAISIE MANUELLE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🎯 SCANNER 3D", "📝 SAISIE MANUELLE"])
    
    with tab1:
        # Simulation simple sans base de données
        st.markdown("### 🎯 SIMULATION SCANNER 3D")
        
        race_selection = st.selectbox("Sélectionnez une race:", 
                                     list(STANDARDS_RACES.keys()),
                                     format_func=lambda x: STANDARDS_RACES[x]['nom_complet'])
        
        if race_selection:
            # Créer des données simulées
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
                        
                        # Afficher quelques points
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
# SECTION 12: PAGE GESTION
# ============================================================================
def page_gestion():
    """Page de gestion du troupeau"""
    st.markdown('<h2 class="section-header">📊 GESTION DU TROUPEAU</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🐑 LISTE", "📈 STATISTIQUES", "🔍 RECHERCHE", "📤 EXPORT"])
    
    with tab1:
        # Données de démonstration
        df = pd.DataFrame({
            'ID': ['HAM-F-2023-001', 'OUDA-F-2023-002', 'SIDAHOU-M-2023-003', 'BERBERE-F-2023-004'],
            'Nom': ['FHAM001', 'FOUDA002', 'MSIDAHOU003', 'FBER004'],
            'Race': ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE'],
            'Sexe': ['F', 'F', 'M', 'F'],
            'Âge': [24, 36, 18, 42],
            'Poids': [55.5, 62.3, 78.9, 48.7],
            'Score': [3, 4, 3, 4],
            'Couleur': ['Rousse', 'Blanche', 'Noire et blanche', 'Brune'],
            'Statut': ['active', 'active', 'active', 'active']
        })
        
        st.dataframe(df, use_container_width=True, height=400)
        
        st.metric("Brebis affichées", len(df))
    
    with tab2:
        st.markdown("### 📊 STATISTIQUES DESCRIPTIVES")
        
        # Statistiques simples
        stats_data = {
            'Race': ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE'],
            'Nombre': [15, 12, 8, 10],
            'Poids moyen (kg)': [58.2, 65.4, 72.1, 45.8],
            'Âge moyen (mois)': [28, 32, 24, 36]
        }
        
        df_stats = pd.DataFrame(stats_data)
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig = px.bar(df_stats, x='Race', y='Poids moyen (kg)',
                        title="Poids moyen par race",
                        color='Nombre',
                        color_continuous_scale='Reds')
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            fig = px.pie(df_stats, values='Nombre', names='Race',
                        title="Répartition des races",
                        hole=0.4)
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 🔍 RECHERCHE SIMPLE")
        
        recherche = st.text_input("Rechercher par nom ou ID")
        
        if recherche:
            st.info(f"Résultats pour: {recherche}")
            # Simulation de résultats
            resultats = pd.DataFrame({
                'ID': [f"{recherche.upper()}-001", f"{recherche.upper()}-002"],
                'Nom': [f"F{recherche.upper()}001", f"M{recherche.upper()}002"],
                'Race': ['HAMRA', 'OUDA'],
                'Poids': [55.5, 62.3]
            })
            st.dataframe(resultats)
    
    with tab4:
        st.markdown("### 📤 EXPORT DE DÉMONSTRATION")
        
        # Créer des données d'exemple
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

# ============================================================================
# SECTION 13: PAGE PRODUCTION
# ============================================================================
def page_production():
    """Page de suivi de production"""
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
                st.success(f"✅ Production enregistrée pour {brebis_id}")
                st.info(f"""
                **Résumé:**
                - Quantité: {quantite} L
                - MG: {mg}%
                - Protéine: {proteine}%
                - Lactose: {lactose}%
                """)
    
    with tab2:
        st.markdown("### 📈 ANALYSE DE PRODUCTION")
        
        # Données simulées
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

# ============================================================================
# SECTION 14: PAGE CRITÈRES DE SÉLECTION
# ============================================================================
def page_criteres():
    """Page des critères de sélection morphologiques et phénotypiques"""
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
            
            # Classification
            st.markdown("### 📊 CLASSIFICATION")
            if score_total >= 4:
                st.success("**Type A (4-5): Excellent pour la production**")
            elif score_total >= 3:
                st.info("**Type B (3-4): Bon pour la production**")
            elif score_total >= 2:
                st.warning("**Type C (2-3): Moyen, à surveiller**")
            else:
                st.error("**Type D (1-2): À améliorer ou réformer**")

# ============================================================================
# SECTION 15: PAGE STATISTIQUES
# ============================================================================
def page_stats():
    """Page d'analyse statistique avancée"""
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

# ============================================================================
# SECTION 16: PAGE GÉNÉTIQUE
# ============================================================================
def page_genetique():
    """Page d'analyse génétique avancée"""
    st.markdown('<h2 class="section-header">🧬 ANALYSE GÉNÉTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🧪 GÉNOTYPAGE", "🌳 DIVERSITÉ"])
    
    with tab1:
        st.markdown("### 🧪 SIMULATION DE GÉNOTYPAGE")
        
        if st.button("🧬 Générer un génotype simulé", type="primary"):
            with st.spinner("Génération en cours..."):
                # Simuler un génotype
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
        
        # Indicateurs simulés
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Ho", "0.325")
        
        with col2:
            st.metric("He", "0.342")
        
        with col3:
            st.metric("Fis", "0.050")

# ============================================================================
# SECTION 17: BARRE LATÉRALE - MODIFIÉE POUR AJOUTER L'OPTION PHOTO
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
         "📸 PHOTO & MESURES",  # ← OPTION AJOUTÉE ICI
         "📐 SCANNER 3D", 
         "📊 GESTION", 
         "🥛 PRODUCTION",
         "🎯 CRITÈRES",
         "📊 RSTATS",
         "🧬 GÉNÉTIQUE"]
    )
    
    st.markdown("---")
    
    # Statistiques rapides
    st.markdown("### 📊 STATISTIQUES")
    st.metric("🐑 Brebis", "20")
    st.metric("🏷️ Races", "6")
    st.metric("🧬 Génotypes", "0")
    
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
        </div>
        """, unsafe_allow_html=True)

# ============================================================================
# SECTION 18: NAVIGATION PRINCIPALE - MODIFIÉE POUR AJOUTER LA PAGE PHOTO
# ============================================================================
if page == "🏠 ACCUEIL":
    page_accueil()
elif page == "📸 PHOTO & MESURES":  # ← PAGE AJOUTÉE ICI
    page_photo_mesures()
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
# SECTION 19: PIED DE PAGE
# ============================================================================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version 5.0</p>
    <p>📸 Photo & Mesures • 📐 Scanner 3D • 🎯 Critères de sélection • 🧬 Génétique • 📊 Statistiques</p>
    <p>© 2024 - Système de gestion scientifique des races ovines algériennes</p>
</div>
""", unsafe_allow_html=True)
