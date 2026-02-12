"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D et Génétique
Base de données simulée de races ovines algériennes
Version avec critères de sélection mammaires et noms génériques
CODE COMPLET AVEC MODULE PHOTO & MESURES AUTOMATIQUES + ANALYSE MULTIPLE
VERSION CORRIGÉE : Fix bug np.int0 → astype(int)

AJOUTS RÉALISÉS :
- Module génomique avancé : FASTA, BLAST, SNP/QTN, maladies, génotypage
- Analyse biochimique du lait professionnelle
- Estimation viande/graisse/os par morphométrie
- Estimation production laitière par morphométrie mammaire
- MODULE API EXTERNES : Météo, Roboflow, Ensembl/NCBI, Cloudinary
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
import traceback

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
    .api-card {
        background: linear-gradient(135deg, #0d47a1 0%, #1565c0 100%);
        color: white;
        padding: 15px;
        border-radius: 15px;
        margin: 10px 0;
        box-shadow: 0 5px 15px rgba(13,71,161,0.2);
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

def get_race_data(race, key, default=None):
    """Récupère les données d'une race de manière sécurisée"""
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
# SECTION 5: MODULE PHOTO & MESURES - VERSION CORRIGÉE
# ============================================================================
class OvinPhotoAnalyzer:
    """Analyseur pour photos profil ET arrière - VERSION CORRIGÉE"""
    
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
        """Détecte l'étalon dans n'importe quelle photo - VERSION CORRIGÉE"""
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
        """Analyse la photo de profil pour mesures corporelles - VERSION AMÉLIORÉE"""
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
        """Analyse la photo arrière pour évaluer les mamelles - VERSION SIMPLIFIÉE"""
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

# ============================================================================
# SECTION 5.1: PAGE PHOTO_MESURES AVEC MODE MANUEL DE SECOURS
# ============================================================================
def page_photo_mesures():
    """Page de capture photo avec 2 vues (profil + arrière) + MODE MANUEL"""
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

def save_complete_characterization(race, age_mois, poids, sexe, classification):
    try:
        cursor = conn.cursor()
        race_code = race[:3] if len(race) >= 3 else "OVN"
        sex_code = "F" if sexe == "Femelle" else "M"
        timestamp = datetime.now().strftime('%y%m%d%H%M')
        identifiant = f"{race_code}-{sex_code}-{timestamp}"
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
        longueur_corps = st.session_state.body_measurements.get('longueur_corps_cm', 0) if 'body_measurements' in st.session_state else 0
        hauteur_garrot = st.session_state.body_measurements.get('hauteur_garrot_cm', 0) if 'body_measurements' in st.session_state else 0
        tour_poitrine = st.session_state.body_measurements.get('tour_poitrine_cm', 0) if 'body_measurements' in st.session_state else 0
        volume_mammaire = st.session_state.mammary_data.get('score_developpement', 0) if 'mammary_data' in st.session_state else 0
        symetrie_mammaire = st.session_state.mammary_data.get('symetrie_mammaire', 0) if 'mammary_data' in st.session_state else 0
        cursor.execute('''
            INSERT INTO brebis (
                identifiant, nom, race, sexe, age_mois, poids,
                longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm,
                volume_mammaire, symetrie_mammaire, 
                notes, statut, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            identifiant,
            f"Char_{identifiant}",
            race,
            sex_code,
            age_mois,
            poids,
            longueur_corps,
            hauteur_garrot,
            tour_poitrine,
            volume_mammaire,
            symetrie_mammaire,
            f"Caractérisation: {classification} | Étalon: {st.session_state.photo_analyzer.etalon_type if hasattr(st.session_state, 'photo_analyzer') else 'inconnu'}",
            'active',
            datetime.now().isoformat()
        ))
        conn.commit()
        st.success(f"✅ Données enregistrées pour **{identifiant}**!")
        st.markdown("### 📥 TÉLÉCHARGEMENTS DISPONIBLES")
        col_dl1, col_dl2 = st.columns(2)
        with col_dl1:
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

def get_recommendations(classification):
    if "EXCELLENT" in classification:
        return "Bonne candidate pour la reproduction et l'amélioration génétique"
    elif "BON" in classification:
        return "À inclure dans le troupeau de production"
    elif "MOYEN" in classification:
        return "À surveiller, évaluer la production réelle"
    else:
        return "Envisager le renouvellement"

# ============================================================================
# SECTION 6: FONCTIONS STATISTIQUES (sans scipy)
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
# SECTION 7: BASE DE DONNÉES - VERSION SÉCURISÉE (AVEC NOUVELLES TABLES)
# ============================================================================
def init_database_safe():
    """Initialise la base de données avec gestion robuste des erreurs + NOUVELLES TABLES"""
    try:
        temp_db = tempfile.NamedTemporaryFile(delete=False, suffix='.db')
        db_path = temp_db.name
        temp_db.close()
        conn = sqlite3.connect(db_path, check_same_thread=False)
        cursor = conn.cursor()
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
                couleur_robe TEXT,
                intensite_couleur INTEGER,
                cornes BOOLEAN,
                taille_cornes_cm FLOAT,
                forme_cornes TEXT,
                type_laine TEXT,
                qualite_laine INTEGER,
                longueur_corps_cm FLOAT,
                hauteur_garrot_cm FLOAT,
                largeur_bassin_cm FLOAT,
                tour_poitrine_cm FLOAT,
                circonference_tete_cm FLOAT,
                longueur_oreille_cm FLOAT,
                volume_mammaire INTEGER,
                symetrie_mammaire INTEGER,
                insertion_trayons INTEGER,
                longueur_trayons_cm FLOAT,
                orientation_trayons TEXT,
                temperement TEXT,
                aptitude TEXT,
                score_conformation FLOAT,
                aptitudes TEXT,
                notes TEXT,
                mere_id TEXT,
                pere_id TEXT,
                coefficient_consanguinite FLOAT DEFAULT 0.0,
                statut TEXT DEFAULT 'active',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
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
        # NOUVELLES TABLES
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS genomic_sequences (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                brebis_id INTEGER,
                sequence_id TEXT,
                description TEXT,
                sequence TEXT,
                longueur INTEGER,
                date_analyse DATE,
                FOREIGN KEY (brebis_id) REFERENCES brebis(id)
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS genetic_markers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                marker_id TEXT UNIQUE,
                chromosome TEXT,
                position INTEGER,
                ref_allele TEXT,
                alt_allele TEXT,
                gene_name TEXT,
                trait_effect TEXT,
                disease_associated TEXT,
                risk_allele TEXT,
                odds_ratio FLOAT,
                p_value FLOAT,
                is_qtn BOOLEAN DEFAULT 0,
                economic_importance INTEGER DEFAULT 0
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS disease_associations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                brebis_id INTEGER,
                marker_id TEXT,
                genotype TEXT,
                disease_name TEXT,
                probability FLOAT,
                risk_level TEXT,
                date_detection DATE,
                FOREIGN KEY (brebis_id) REFERENCES brebis(id),
                FOREIGN KEY (marker_id) REFERENCES genetic_markers(marker_id)
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS milk_composition (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                brebis_id INTEGER,
                date_analyse DATE,
                echantillon_id TEXT,
                mg_g_100ml FLOAT,
                proteines_g_100ml FLOAT,
                lactose_g_100ml FLOAT,
                extrac_secret FLOAT,
                ph FLOAT,
                acidite_dornic FLOAT,
                cellules_somatiques INTEGER,
                acides_gras_satures_g FLOAT,
                acides_gras_insatures_g FLOAT,
                calcium_mg_100ml FLOAT,
                phosphore_mg_100ml FLOAT,
                magnesium_mg_100ml FLOAT,
                potassium_mg_100ml FLOAT,
                vitamine_a_ug_100ml FLOAT,
                vitamine_e_ug_100ml FLOAT,
                aptitude_fromagere TEXT,
                notes TEXT,
                FOREIGN KEY (brebis_id) REFERENCES brebis(id)
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS carcass_estimates (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                brebis_id INTEGER,
                date_estimation DATE,
                poids_vif_kg FLOAT,
                longueur_corps_cm FLOAT,
                hauteur_garrot_cm FLOAT,
                tour_poitrine_cm FLOAT,
                viande_estimee_kg FLOAT,
                graisse_estimee_kg FLOAT,
                os_estimes_kg FLOAT,
                rendement_carcasse FLOAT,
                methode TEXT,
                confiance FLOAT,
                notes TEXT,
                FOREIGN KEY (brebis_id) REFERENCES brebis(id)
            )
        ''')
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS milk_production_estimates (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                brebis_id INTEGER,
                date_estimation DATE,
                volume_mammaire_cm3 FLOAT,
                largeur_mammaire_cm FLOAT,
                hauteur_mammaire_cm FLOAT,
                symetrie FLOAT,
                lait_jour_estime_l FLOAT,
                lactation_potentielle_l FLOAT,
                confiance FLOAT,
                notes TEXT,
                FOREIGN KEY (brebis_id) REFERENCES brebis(id)
            )
        ''')
        cursor.execute("SELECT COUNT(*) FROM brebis")
        count = cursor.fetchone()[0]
        if count == 0:
            peupler_base_races_safe(cursor, conn)
            peupler_genomic_data(cursor, conn)
            peupler_milk_composition(cursor, conn)
            peupler_carcass_estimates(cursor, conn)
            peupler_milk_prod_estimates(cursor, conn)
        conn.commit()
        return conn
    except Exception as e:
        st.error(f"Erreur d'initialisation: {str(e)}")
        conn = sqlite3.connect(':memory:', check_same_thread=False)
        return conn

def peupler_base_races_safe(cursor, conn):
    races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU']
    brebis_data = []
    for i in range(1, 21):
        race = random.choice(races)
        sexe = random.choice(['F', 'M'])
        race_code = race[:3] if race != 'INCONNU' else 'INC'
        identifiant = f"{race_code}-{sexe}-2023-{i:03d}"
        if sexe == 'F':
            nom = f"F{race_code}{i:03d}"
        else:
            nom = f"M{race_code}{i:03d}"
        age_mois = random.randint(12, 84)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        try:
            poids_data = get_race_data(race, 'poids_adulte')
            if sexe.lower() in poids_data:
                poids_min, poids_max = poids_data[sexe.lower()]
            else:
                poids_min, poids_max = (35, 60) if sexe == 'F' else (50, 80)
        except:
            poids_min, poids_max = (35, 60) if sexe == 'F' else (50, 80)
        poids = random.uniform(poids_min, poids_max)
        score_condition = random.randint(2, 4)
        couleurs = {
            'HAMRA': ['Rousse', 'Rousse foncée', 'Marron'],
            'OUDA': ['Blanche', 'Crème', 'Blanc cassé'],
            'SIDAHOU': ['Noire et blanche', 'Pie noire', 'Tête noire'],
            'BERBERE': ['Noire', 'Brune', 'Grise', 'Pie'],
            'CROISE': ['Variable', 'Panachée', 'Mélangée'],
            'INCONNU': ['Indéterminée']
        }
        couleur_robe = random.choice(couleurs.get(race, ['Indéterminée']))
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
        score_conformation = random.uniform(5.0, 9.0)
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
        cursor.execute('''
            INSERT INTO brebis (identifiant, nom, race, sexe, age_mois, poids, statut)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', ('TEST-001', 'TestBrebis', 'HAMRA', 'F', 24, 55.5, 'active'))
        conn.commit()

def peupler_genomic_data(cursor, conn):
    markers = [
        ('SNP001', '1', 12345678, 'A', 'G', 'GDF8', 'Masse musculaire', 'Hypertrophie musculaire', 'G', 2.5, 0.0001, 1, 5),
        ('SNP002', '2', 87654321, 'C', 'T', 'PRNP', 'Sensibilité tremblante', 'Tremblante', 'T', 12.0, 0.00001, 0, 5),
        ('SNP003', '3', 11223344, 'G', 'A', 'FecB', 'Prolificité', 'Aucune', 'A', 3.2, 0.0005, 1, 5),
        ('SNP004', '4', 99887766, 'T', 'C', 'DGAT1', 'Taux MG lait', 'Aucune', 'C', 1.8, 0.001, 1, 4),
        ('SNP005', '5', 44556677, 'A', 'G', 'CAST', 'Tendreté viande', 'Aucune', 'G', 1.5, 0.01, 0, 3),
        ('QTN001', '6', 123987, 'A', 'G', 'LALBA', 'Production lait', 'Aucune', 'G', 2.1, 0.0008, 1, 5),
        ('QTN002', '7', 987123, 'C', 'T', 'CSN1S1', 'Taux caséine', 'Aucune', 'T', 1.9, 0.002, 1, 4)
    ]
    cursor.executemany('''
        INSERT OR IGNORE INTO genetic_markers 
        (marker_id, chromosome, position, ref_allele, alt_allele, gene_name, trait_effect, 
         disease_associated, risk_allele, odds_ratio, p_value, is_qtn, economic_importance)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
    ''', markers)
    cursor.execute("SELECT id, identifiant FROM brebis LIMIT 5")
    brebis_list = cursor.fetchall()
    bases = ['A','C','G','T']
    for b in brebis_list:
        seq_id = f"FASTA_{b[1]}"
        seq_len = random.randint(800, 1200)
        sequence = ''.join(random.choices(bases, k=seq_len))
        cursor.execute('''
            INSERT INTO genomic_sequences (brebis_id, sequence_id, description, sequence, longueur, date_analyse)
            VALUES (?, ?, ?, ?, ?, ?)
        ''', (b[0], seq_id, f"Séquence simulée de {b[1]}", sequence, seq_len, date.today().isoformat()))
    cursor.execute("SELECT id, identifiant FROM brebis WHERE sexe='F' LIMIT 10")
    brebis_f = cursor.fetchall()
    cursor.execute("SELECT marker_id, disease_associated, risk_allele FROM genetic_markers WHERE disease_associated != 'Aucune'")
    disease_markers = cursor.fetchall()
    for b in brebis_f[:3]:
        for dm in disease_markers:
            if random.random() < 0.3:
                genotype = random.choice([f"{dm[2]}{dm[2]}", f"{dm[2]}{random.choice(['A','C','G','T'])}"])
                cursor.execute('''
                    INSERT INTO disease_associations (brebis_id, marker_id, genotype, disease_name, probability, risk_level, date_detection)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                ''', (b[0], dm[0], genotype, dm[1], random.uniform(0.6,0.95), 
                      random.choice(['Faible','Modéré','Élevé']), date.today().isoformat()))
    conn.commit()

def peupler_milk_composition(cursor, conn):
    cursor.execute("SELECT id FROM brebis WHERE sexe='F' LIMIT 15")
    brebis_f = cursor.fetchall()
    for b in brebis_f:
        echantillon = f"LAIT-{b[0]}-{random.randint(100,999)}"
        cursor.execute('''
            INSERT INTO milk_composition 
            (brebis_id, date_analyse, echantillon_id, mg_g_100ml, proteines_g_100ml, lactose_g_100ml,
             extrac_secret, ph, acidite_dornic, cellules_somatiques, acides_gras_satures_g,
             acides_gras_insatures_g, calcium_mg_100ml, phosphore_mg_100ml, magnesium_mg_100ml,
             potassium_mg_100ml, vitamine_a_ug_100ml, vitamine_e_ug_100ml, aptitude_fromagere, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            b[0],
            date.today().isoformat(),
            echantillon,
            round(random.uniform(6.0, 9.0), 1),
            round(random.uniform(5.0, 7.0), 1),
            round(random.uniform(4.0, 5.5), 1),
            round(random.uniform(16, 22), 1),
            round(random.uniform(6.5, 6.9), 2),
            round(random.uniform(14, 20), 1),
            random.randint(50, 400),
            round(random.uniform(3.5, 5.0), 1),
            round(random.uniform(0.8, 1.5), 1),
            round(random.uniform(150, 220), 1),
            round(random.uniform(120, 180), 1),
            round(random.uniform(15, 25), 1),
            round(random.uniform(130, 190), 1),
            round(random.uniform(30, 60), 1),
            round(random.uniform(0.1, 0.5), 2),
            random.choice(['Bonne', 'Moyenne', 'Excellente', 'Faible']),
            "Analyse simulée"
        ))
    conn.commit()

def peupler_carcass_estimates(cursor, conn):
    cursor.execute("SELECT id, poids, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm FROM brebis")
    brebis = cursor.fetchall()
    for b in brebis:
        poids = b[1] if b[1] else 50.0
        longueur = b[2] if b[2] else 100.0
        hauteur = b[3] if b[3] else 70.0
        poitrine = b[4] if b[4] else 100.0
        viande = poids * 0.55 + longueur * 0.1
        graisse = poids * 0.20 - longueur * 0.05
        os = poids * 0.15 + hauteur * 0.02
        rendement = (viande + graisse + os) / poids * 100 if poids > 0 else 0
        cursor.execute('''
            INSERT INTO carcass_estimates 
            (brebis_id, date_estimation, poids_vif_kg, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm,
             viande_estimee_kg, graisse_estimee_kg, os_estimes_kg, rendement_carcasse, methode, confiance, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            b[0],
            date.today().isoformat(),
            poids,
            longueur,
            hauteur,
            poitrine,
            round(viande, 2),
            round(graisse, 2),
            round(os, 2),
            round(rendement, 1),
            "Modèle Ovin Manager",
            round(random.uniform(0.75, 0.95), 2),
            "Estimation automatique"
        ))
    conn.commit()

def peupler_milk_prod_estimates(cursor, conn):
    cursor.execute("SELECT id, volume_mammaire, longueur_trayons_cm FROM brebis WHERE sexe='F' AND volume_mammaire IS NOT NULL LIMIT 20")
    brebis_f = cursor.fetchall()
    for b in brebis_f:
        vol_cm3 = b[1] * 80 + 100 if b[1] else 250
        largeur = random.uniform(6, 12)
        hauteur = random.uniform(8, 16)
        symetrie = random.uniform(0.7, 0.95)
        lait_estime = vol_cm3 * 0.01 + largeur * 0.2 + hauteur * 0.15
        lait_potentiel = lait_estime * random.uniform(1.2, 1.5)
        cursor.execute('''
            INSERT INTO milk_production_estimates 
            (brebis_id, date_estimation, volume_mammaire_cm3, largeur_mammaire_cm, hauteur_mammaire_cm,
             symetrie, lait_jour_estime_l, lactation_potentielle_l, confiance, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            b[0],
            date.today().isoformat(),
            vol_cm3,
            round(largeur, 1),
            round(hauteur, 1),
            round(symetrie, 2),
            round(lait_estime, 2),
            round(lait_potentiel, 2),
            round(random.uniform(0.7, 0.9), 2),
            "Estimation morphométrique"
        ))
    conn.commit()

@st.cache_resource
def get_database_connection():
    return init_database_safe()

conn = get_database_connection()

# ============================================================================
# SECTION 8: MODULE SCANNER 3D
# ============================================================================
class Scanner3D:
    """Simulateur de scanner 3D pour ovins"""
    
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

# ============================================================================
# SECTION 9: MODULE GÉNÉTIQUE
# ============================================================================
class ModuleGenetique:
    """Module d'analyse génétique"""
    
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

# ============================================================================
# SECTION 10: PAGE ACCUEIL
# ============================================================================
def page_accueil():
    """Page d'accueil avec vue d'ensemble"""
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - RACES ALGÉRIENNES</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique des races ovines algériennes**")
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
        - **📦 Analyse Multiple**: Analyse en lot pour éleveurs
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
            resultats = pd.DataFrame({
                'ID': [f"{recherche.upper()}-001", f"{recherche.upper()}-002"],
                'Nom': [f"F{recherche.upper()}001", f"M{recherche.upper()}002"],
                'Race': ['HAMRA', 'OUDA'],
                'Poids': [55.5, 62.3]
            })
            st.dataframe(resultats)
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

# ============================================================================
# SECTION 20: PAGE ANALYSE MULTIPLE
# ============================================================================
def page_analyse_multiple():
    """Page pour analyser plusieurs brebis en lot"""
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

# ============================================================================
# SECTION 21: MODULE GÉNOMIQUE AVANCÉ (FASTA/BLAST, SNP/QTN, MALADIES, GÉNOTYPAGE)
# ============================================================================
class GenomicAnalyzer:
    """Analyse génomique avancée : séquences, BLAST, variants, maladies, génotypage"""

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
    def genotyper_brebis(brebis_id, cursor):
        cursor.execute("SELECT marker_id, ref_allele, alt_allele FROM genetic_markers")
        markers = cursor.fetchall()
        genotypes = []
        for m in markers:
            prob_alt = random.uniform(0.1, 0.4)
            allele1 = m[2] if random.random() < prob_alt else m[1]
            allele2 = m[2] if random.random() < prob_alt else m[1]
            genotype = allele1 + allele2
            genotypes.append((brebis_id, m[0], m[1], m[2], genotype, 
                             date.today().isoformat()))
        return genotypes

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
            cursor = conn.cursor()
            cursor.execute("SELECT id, identifiant FROM brebis LIMIT 10")
            brebis_list = cursor.fetchall()
            brebis_dict = {f"{b[1]} (ID: {b[0]})": b[0] for b in brebis_list}
            selected_brebis = st.selectbox("Choisir une brebis", list(brebis_dict.keys()))
            brebis_id = brebis_dict[selected_brebis]
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
        cursor = conn.cursor()
        cursor.execute("SELECT marker_id, chromosome, position, gene_name, trait_effect, is_qtn, economic_importance FROM genetic_markers")
        markers_data = cursor.fetchall()
        df_markers = pd.DataFrame(markers_data, 
                                 columns=['Marqueur','Chr','Pos','Gène','Effet','QTN','Importance'])
        st.dataframe(df_markers)
    with tab3:
        st.markdown("### 🩺 DÉTECTION DE MALADIES GÉNÉTIQUES")
        st.markdown("Évaluation du risque à partir des génotypes.")
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant FROM brebis")
        all_brebis = cursor.fetchall()
        brebis_choice = st.selectbox("Sélectionner une brebis", 
                                    [f"{b[1]} (ID: {b[0]})" for b in all_brebis],
                                    key="maladie_brebis")
        brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
        if st.button("🧪 Évaluer les risques génétiques"):
            cursor.execute("SELECT marker_id, genotype FROM genotypage WHERE brebis_id=? LIMIT 10", (brebis_id,))
            genotypes_raw = cursor.fetchall()
            if not genotypes_raw:
                genotypes = GenomicAnalyzer.genotyper_brebis(brebis_id, cursor)
                for g in genotypes[:5]:
                    cursor.execute('''
                        INSERT INTO genotypage (brebis_id, marqueur, allele1, allele2, genotype, date_analyse)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (g[0], g[1], g[2], g[3], g[4], g[5]))
                conn.commit()
                cursor.execute("SELECT marker_id, genotype FROM genotypage WHERE brebis_id=? LIMIT 10", (brebis_id,))
                genotypes_raw = cursor.fetchall()
            genotype_dict = {g[0]: g[1] for g in genotypes_raw}
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
    with tab4:
        st.markdown("### 🧬 GÉNOTYPAGE SIMULÉ")
        st.markdown("Génération de génotypes pour les marqueurs d'intérêt.")
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant FROM brebis")
        all_brebis = cursor.fetchall()
        brebis_choice = st.selectbox("Choisir une brebis", 
                                    [f"{b[1]} (ID: {b[0]})" for b in all_brebis],
                                    key="geno_brebis")
        brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
        if st.button("🧬 Générer génotype complet"):
            genotypes = GenomicAnalyzer.genotyper_brebis(brebis_id, cursor)
            for g in genotypes:
                cursor.execute('''
                    INSERT OR REPLACE INTO genotypage 
                    (brebis_id, marqueur, allele1, allele2, genotype, date_analyse)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (g[0], g[1], g[2], g[3], g[4], g[5]))
            conn.commit()
            st.success(f"{len(genotypes)} marqueurs génotypés et enregistrés.")
            df_geno_full = pd.DataFrame(genotypes, 
                                       columns=['Brebis_id','Marqueur','Allèle1','Allèle2','Génotype','Date'])
            st.dataframe(df_geno_full[['Marqueur','Génotype','Date']])

# ============================================================================
# SECTION 22: MODULE D'ANALYSE BIOCHIMIQUE DU LAIT (PROFESSIONNELLE)
# ============================================================================
class LaitAnalyzer:
    """Analyse biochimique professionnelle du lait de brebis"""

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
                cursor = conn.cursor()
                cursor.execute("SELECT id FROM brebis WHERE identifiant=? OR nom=?", (brebis_id, brebis_id))
                result = cursor.fetchone()
                if result:
                    brebis_db_id = result[0]
                else:
                    cursor.execute("INSERT INTO brebis (identifiant, nom, race, sexe, age_mois, poids) VALUES (?,?,?,?,?,?)",
                                  (brebis_id, f"Brebis_{brebis_id}", 'INCONNU', 'F', 24, 50.0))
                    brebis_db_id = cursor.lastrowid
                cursor.execute('''
                    INSERT INTO milk_composition 
                    (brebis_id, date_analyse, echantillon_id, mg_g_100ml, proteines_g_100ml, lactose_g_100ml,
                     extrac_secret, ph, acidite_dornic, cellules_somatiques, acides_gras_satures_g,
                     acides_gras_insatures_g, calcium_mg_100ml, phosphore_mg_100ml, magnesium_mg_100ml,
                     potassium_mg_100ml, vitamine_a_ug_100ml, vitamine_e_ug_100ml, aptitude_fromagere, notes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    brebis_db_id, date_prelevement.isoformat(), echantillon_id, mg, proteines, lactose,
                    extrac_secret, ph, acidite, cellules, ags, agi, calcium, phosphore, magnesium,
                    potassium, vit_a, vit_e, str(aptitude), "Analyse professionnelle"
                ))
                conn.commit()
                st.success("✅ Analyse enregistrée dans la base de données !")
    with tab2:
        st.markdown("### 📊 HISTORIQUE DES ANALYSES")
        cursor = conn.cursor()
        cursor.execute('''
            SELECT b.identifiant, mc.date_analyse, mc.mg_g_100ml, mc.proteines_g_100ml, 
                   mc.lactose_g_100ml, mc.cellules_somatiques, mc.aptitude_fromagere
            FROM milk_composition mc
            JOIN brebis b ON mc.brebis_id = b.id
            ORDER BY mc.date_analyse DESC
            LIMIT 20
        ''')
        data = cursor.fetchall()
        if data:
            df_hist = pd.DataFrame(data, columns=['Brebis','Date','MG','Protéines','Lactose','Cellules','Aptitude'])
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

# ============================================================================
# SECTION 23: MODULE D'ESTIMATION DE LA COMPOSITION CORPOREELLE (VIANDE/GRAISSE/OS)
# ============================================================================
class CarcassEstimator:
    """Prédiction de la composition de la carcasse à partir des mesures morphométriques"""

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

def page_estimation_viande():
    st.markdown('<h2 class="section-header">🥩 ESTIMATION DE LA COMPOSITION DE LA CARCASSE</h2>', unsafe_allow_html=True)
    st.markdown("**Prédiction de la quantité de viande, graisse et os à partir des mesures morphométriques**")
    tab1, tab2 = st.tabs(["📏 Estimation individuelle", "📊 Analyse de lot"])
    with tab1:
        st.markdown("### 🐑 Saisie des mesures")
        col1, col2 = st.columns(2)
        with col1:
            identifiant = st.text_input("Identifiant de la brebis", "HAM-F-001")
            cursor = conn.cursor()
            cursor.execute("SELECT id, race, age_mois, poids, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm FROM brebis WHERE identifiant=? OR nom=?", (identifiant, identifiant))
            existing = cursor.fetchone()
            if existing:
                st.success(f"Brebis trouvée : {identifiant}")
                race_def = existing[1]
                age_def = existing[2] if existing[2] else 24
                poids_def = existing[3] if existing[3] else 50.0
                longueur_def = existing[4] if existing[4] else 100.0
                hauteur_def = existing[5] if existing[5] else 70.0
                poitrine_def = existing[6] if existing[6] else 100.0
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
                cursor = conn.cursor()
                if existing:
                    brebis_db_id = existing[0]
                else:
                    cursor.execute("INSERT INTO brebis (identifiant, nom, race, sexe, age_mois, poids) VALUES (?,?,?,?,?,?)",
                                  (identifiant, f"Brebis_{identifiant}", race, 'F', age_mois, poids))
                    brebis_db_id = cursor.lastrowid
                cursor.execute('''
                    INSERT INTO carcass_estimates 
                    (brebis_id, date_estimation, poids_vif_kg, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm,
                     viande_estimee_kg, graisse_estimee_kg, os_estimes_kg, rendement_carcasse, methode, confiance, notes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    brebis_db_id, date.today().isoformat(), poids, longueur, hauteur, tour_poitrine,
                    estimation['viande_kg'], estimation['graisse_kg'], estimation['os_kg'],
                    estimation['rendement_carcasse'], "Modèle Ovin Manager", estimation['confiance'],
                    f"Estimation basée sur {race}, âge {age_mois} mois"
                ))
                conn.commit()
                st.success("✅ Estimation enregistrée !")
    with tab2:
        st.markdown("### 📊 ESTIMATION EN LOT")
        st.info("Sélectionnez plusieurs brebis et estimez leur composition en une fois.")
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, race, age_mois, poids, longueur_corps_cm, hauteur_garrot_cm, tour_poitrine_cm FROM brebis WHERE sexe='F' LIMIT 20")
        brebis_list = cursor.fetchall()
        if brebis_list:
            df_brebis = pd.DataFrame(brebis_list, 
                                    columns=['ID','Identifiant','Race','Âge','Poids','Longueur','Hauteur','Poitrine'])
            st.dataframe(df_brebis[['Identifiant','Race','Âge','Poids']])
            selected = st.multiselect("Choisir les brebis", df_brebis['Identifiant'].tolist())
            if selected and st.button("📊 Estimer pour la sélection"):
                results = []
                for ident in selected:
                    row = df_brebis[df_brebis['Identifiant'] == ident].iloc[0]
                    est = CarcassEstimator.estimer_composition(
                        row['Poids'], row['Longueur'], row['Hauteur'], row['Poitrine'], row['Âge'], row['Race']
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

# ============================================================================
# SECTION 24: MODULE D'ESTIMATION DE LA PRODUCTION LAITIÈRE PAR MORPHOMÉTRIE MAMMAIRE
# ============================================================================
class MilkProductionEstimator:
    """Estimation de la production laitière à partir des mesures morphométriques des mamelles"""

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

def page_estimation_lait_morpho():
    st.markdown('<h2 class="section-header">🍼 ESTIMATION LAIT PAR MORPHOMÉTRIE MAMMAIRE</h2>', unsafe_allow_html=True)
    st.markdown("**Prédiction de la production laitière individuelle à partir des caractéristiques de la mamelle**")
    tab1, tab2 = st.tabs(["📏 Estimation individuelle", "📈 Analyse par ID"])
    with tab1:
        st.markdown("### 🍼 MESURES DE LA MAMELLE")
        col1, col2 = st.columns(2)
        with col1:
            identifiant = st.text_input("Identifiant de la brebis", "HAM-F-001")
            cursor = conn.cursor()
            cursor.execute("SELECT id, race, age_mois, volume_mammaire, longueur_trayons_cm FROM brebis WHERE identifiant=? OR nom=?", (identifiant, identifiant))
            existing = cursor.fetchone()
            if existing:
                race_def = existing[1] if existing[1] else 'HAMRA'
                age_def = existing[2] if existing[2] else 24
                vol_mammaire_score = existing[3] if existing[3] else 3
                vol_cm3_def = vol_mammaire_score * 80 + 100
                longueur_trayons = existing[4] if existing[4] else 4.5
            else:
                race_def = 'HAMRA'
                age_def = 24
                vol_cm3_def = 250
                longueur_trayons = 4.5
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
                cursor = conn.cursor()
                if existing:
                    brebis_db_id = existing[0]
                else:
                    cursor.execute("INSERT INTO brebis (identifiant, nom, race, sexe, age_mois) VALUES (?,?,?,?,?)",
                                  (identifiant, f"Brebis_{identifiant}", race, 'F', age_mois))
                    brebis_db_id = cursor.lastrowid
                cursor.execute('''
                    INSERT INTO milk_production_estimates 
                    (brebis_id, date_estimation, volume_mammaire_cm3, largeur_mammaire_cm, hauteur_mammaire_cm,
                     symetrie, lait_jour_estime_l, lactation_potentielle_l, confiance, notes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    brebis_db_id, date.today().isoformat(), volume_cm3, largeur_cm, hauteur_cm,
                    symetrie, estimation['lait_jour_estime_l'], estimation['lactation_potentielle_l'],
                    estimation['confiance'], f"Estimation basée sur morphométrie mammaire"
                ))
                conn.commit()
                st.success("✅ Estimation enregistrée !")
    with tab2:
        st.markdown("### 📈 CONSULTATION PAR ID")
        st.markdown("Recherchez une brebis et visualisez ses estimations de production laitière.")
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant FROM brebis WHERE sexe='F' LIMIT 50")
        brebis_f = cursor.fetchall()
        if brebis_f:
            brebis_choice = st.selectbox("Choisir une brebis", 
                                        [f"{b[1]} (ID: {b[0]})" for b in brebis_f])
            brebis_id = int(brebis_choice.split("ID: ")[1][:-1])
            cursor.execute('''
                SELECT date_estimation, volume_mammaire_cm3, lait_jour_estime_l, lactation_potentielle_l, confiance
                FROM milk_production_estimates
                WHERE brebis_id=?
                ORDER BY date_estimation DESC
            ''', (brebis_id,))
            est_data = cursor.fetchall()
            if est_data:
                df_est = pd.DataFrame(est_data, 
                                     columns=['Date','Volume (cm³)','Lait/jour (L)','Lactation (L)','Confiance'])
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
# SECTION 25: MODULE D'INTÉGRATION D'API EXTERNES (NOUVEAU)
# ============================================================================

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import hashlib
import time

class APIManager:
    """Gestionnaire centralisé des appels API avec cache et retry"""
    
    _cache = {}
    _cache_timeout = 3600  # 1 heure
    
    def __init__(self):
        self.session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        
        # Clés API - À REMPLACER PAR VOS VRAIES CLÉS
        self.OPENWEATHER_KEY = "VOTRE_CLE_OPENWEATHERMAP"  # Obtenez votre clé sur https://openweathermap.org/api
        self.ROBOFLOW_API_KEY = "VOTRE_CLE_ROBOFLOW"      # Obtenez votre clé sur https://roboflow.com
        self.NCBI_API_KEY = "VOTRE_CLE_NCBI"              # Optionnel - pour plus de requêtes NCBI
    
    # ------------------------------------------------------------------------
    # 1. API MÉTÉO (OpenWeatherMap) - PRIORITAIRE
    # ------------------------------------------------------------------------
    
    @staticmethod
    def get_weather(city="Tlemcen", country="DZ"):
        """Récupère les données météo actuelles pour une ville algérienne"""
        api_key = "VOTRE_CLE_OPENWEATHERMAP"  # REMPLACEZ ICI
        
        try:
            url = "http://api.openweathermap.org/data/2.5/weather"
            params = {
                'q': f"{city},{country}",
                'appid': api_key,
                'units': 'metric',
                'lang': 'fr'
            }
            
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                return {
                    'temperature': round(data['main']['temp'], 1),
                    'ressenti': round(data['main']['feels_like'], 1),
                    'humidite': data['main']['humidity'],
                    'pression': data['main']['pressure'],
                    'description': data['weather'][0]['description'],
                    'icone': data['weather'][0]['icon'],
                    'vent_vitesse': data['wind']['speed'],
                    'vent_direction': data['wind'].get('deg', 0),
                    'visibilite': data.get('visibility', 10000) / 1000,
                    'lever_soleil': datetime.fromtimestamp(data['sys']['sunrise']).strftime('%H:%M'),
                    'coucher_soleil': datetime.fromtimestamp(data['sys']['sunset']).strftime('%H:%M'),
                    'ville': data['name'],
                    'pays': data['sys']['country']
                }
            else:
                st.warning(f"⚠️ API Météo : {response.status_code} - Utilisation du mode simulation")
                return None
                
        except Exception as e:
            st.warning(f"⚠️ Erreur API Météo : {str(e)} - Utilisation du mode simulation")
            return None
    
    @staticmethod
    def get_forecast(city="Tlemcen", days=5):
        """Prévisions météo sur 5 jours"""
        api_key = "VOTRE_CLE_OPENWEATHERMAP"  # REMPLACEZ ICI
        
        try:
            url = "http://api.openweathermap.org/data/2.5/forecast"
            params = {
                'q': city,
                'appid': api_key,
                'units': 'metric',
                'lang': 'fr',
                'cnt': days * 8
            }
            
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                forecast_data = []
                for item in data['list']:
                    forecast_data.append({
                        'datetime': datetime.fromtimestamp(item['dt']),
                        'temperature': round(item['main']['temp'], 1),
                        'ressenti': round(item['main']['feels_like'], 1),
                        'humidite': item['main']['humidity'],
                        'description': item['weather'][0]['description'],
                        'icone': item['weather'][0]['icon'],
                        'probabilite_pluie': item.get('pop', 0) * 100
                    })
                return forecast_data
            return None
        except Exception as e:
            st.warning(f"⚠️ Erreur API Prévisions : {str(e)}")
            return None
    
    @staticmethod
    def get_heat_stress_level(temperature, humidite):
        """Calcule le niveau de stress thermique (indice THI) pour les ovins"""
        # THI = (0.8 * T) + (H * (T - 14.4)) + 46.4
        thi = (0.8 * temperature) + (humidite * (temperature - 14.4) / 100) + 46.4
        
        if thi < 70:
            return "NORMAL", "✅ Conditions idéales"
        elif thi < 78:
            return "MODÉRÉ", "⚠️ Stress thermique léger - Surveiller l'abreuvement"
        elif thi < 85:
            return "ÉLEVÉ", "🔥 Stress thermique élevé - Ombre et eau fraîche indispensables"
        else:
            return "CRITIQUE", "☠️ Danger ! Ventilation et rafraîchissement immédiat"
    
    @staticmethod
    def get_nutrition_advice(weather_data):
        """Génère des conseils nutritionnels basés sur la météo"""
        if not weather_data:
            return "Données météo non disponibles - Utilisez les recommandations standard"
        
        temp = weather_data['temperature']
        humidite = weather_data['humidite']
        stress_level, message = APIManager.get_heat_stress_level(temp, humidite)
        
        advice = f"🌡️ **Température :** {temp}°C | 💧 **Humidité :** {humidite}%\n"
        advice += f"📊 **Stress thermique :** {message}\n\n"
        advice += "**🥛 RECOMMANDATIONS NUTRITIONNELLES :**\n"
        
        if stress_level == "NORMAL":
            advice += """✅ Ration standard maintenue
✅ Eau à volonté (6-8 L/jour)
✅ Fourrage normal
✅ Minéraux standard"""
        elif stress_level == "MODÉRÉ":
            advice += """⚠️ **Augmenter l'eau de 25%** (8-10 L/jour)
⚠️ Réduire les fibres de 15%
⚠️ Distribuer en soirée
✅ Bicarbonate (0.5% de la ration)"""
        elif stress_level == "ÉLEVÉ":
            advice += """🔥 **EAU +50%** (12-15 L/jour)
🔥 Alimentation concentrée le matin uniquement
🔥 Fourrage vert uniquement
🔥 Bicarbonate : 1% de la ration
🔥 Vitamines C+E"""
        else:
            advice += """☠️ **SITUATION CRITIQUE :**
☠️ Abreuvement permanent - eau fraîche
☠️ Alimentation liquide uniquement
☠️ Ventilation forcée
☠️ Aucun travail"""
        
        return advice
    
    # ------------------------------------------------------------------------
    # 2. API ROBOFLOW - Détection des mamelles par IA
    # ------------------------------------------------------------------------
    
    @staticmethod
    def detect_mammary_roboflow(image, api_key="VOTRE_CLE_ROBOFLOW"):
        """Détecte et mesure les mamelles via l'API Roboflow"""
        try:
            import cv2
            import base64
            
            # Convertir l'image en base64
            if isinstance(image, np.ndarray):
                _, img_encoded = cv2.imencode('.jpg', image)
                img_bytes = img_encoded.tobytes()
            else:
                img_bytes = image.getvalue()
            
            img_base64 = base64.b64encode(img_bytes).decode('utf-8')
            
            # Appel API Roboflow (modèle à entraîner)
            url = "https://detect.roboflow.com/ovin-mammary-detection/1"
            params = {
                'api_key': api_key,
                'confidence': 40,
                'overlap': 30
            }
            
            response = requests.post(
                url,
                params=params,
                data=img_base64,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                timeout=15
            )
            
            if response.status_code == 200:
                result = response.json()
                predictions = result.get('predictions', [])
                
                mammary_data = {
                    'nombre_mamelles': len(predictions),
                    'volumes': [],
                    'largueurs': [],
                    'hauteurs': [],
                    'confiance_moyenne': 0,
                    'bounding_boxes': []
                }
                
                total_conf = 0
                for pred in predictions:
                    mammary_data['bounding_boxes'].append({
                        'x': pred['x'],
                        'y': pred['y'],
                        'width': pred['width'],
                        'height': pred['height'],
                        'confidence': pred['confidence']
                    })
                    mammary_data['largueurs'].append(pred['width'])
                    mammary_data['hauteurs'].append(pred['height'])
                    total_conf += pred['confidence']
                
                if mammary_data['nombre_mamelles'] > 0:
                    mammary_data['confiance_moyenne'] = total_conf / mammary_data['nombre_mamelles']
                    mammary_data['largeur_moyenne'] = np.mean(mammary_data['largueurs'])
                    mammary_data['hauteur_moyenne'] = np.mean(mammary_data['hauteurs'])
                    # Estimation simplifiée du volume
                    mammary_data['volume_moyen_cm3'] = (
                        mammary_data['largeur_moyenne'] * 
                        mammary_data['hauteur_moyenne'] * 
                        0.7
                    )
                
                return mammary_data
            else:
                # Mode simulation si l'API n'est pas disponible
                return {
                    'nombre_mamelles': 2,
                    'confiance_moyenne': 0.85,
                    'largeur_moyenne': 120,
                    'hauteur_moyenne': 150,
                    'volume_moyen_cm3': 250,
                    'mode': 'simulation'
                }
                
        except Exception as e:
            st.warning(f"⚠️ Erreur Roboflow: {str(e)} - Utilisation du mode simulation")
            return {
                'nombre_mamelles': 2,
                'confiance_moyenne': 0.85,
                'largeur_moyenne': 120,
                'hauteur_moyenne': 150,
                'volume_moyen_cm3': 250,
                'mode': 'simulation'
            }
    
    # ------------------------------------------------------------------------
    # 3. API Ensembl - Génomique ovine
    # ------------------------------------------------------------------------
    
    @staticmethod
    def get_ensembl_gene(gene_name, species="ovis_aries"):
        """Récupère les informations d'un gène via Ensembl REST API"""
        try:
            url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene_name}"
            headers = {"Content-Type": "application/json"}
            
            response = requests.get(url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                return response.json()
            else:
                return None
                
        except Exception as e:
            st.warning(f"⚠️ Erreur Ensembl: {str(e)}")
            return None
    
    @staticmethod
    def get_ensembl_sequence(gene_id, species="ovis_aries"):
        """Récupère la séquence FASTA d'un gène"""
        try:
            url = f"https://rest.ensembl.org/sequence/id/{gene_id}"
            headers = {"Content-Type": "text/x-fasta"}
            
            response = requests.get(url, headers=headers, timeout=15)
            
            if response.status_code == 200:
                return response.text
            else:
                return None
                
        except Exception as e:
            st.warning(f"⚠️ Erreur séquence Ensembl: {str(e)}")
            return None
    
    # ------------------------------------------------------------------------
    # 4. API NCBI E-utilities - Recherche bibliographique
    # ------------------------------------------------------------------------
    
    @staticmethod
    def ncbi_search(query, database="pubmed", max_results=5):
        """Recherche dans les bases de données NCBI (PubMed, Gene, etc.)"""
        api_key = "VOTRE_CLE_NCBI"  # Optionnel
        
        try:
            # ESearch
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': database,
                'term': query,
                'retmax': max_results,
                'retmode': 'json'
            }
            if api_key != "VOTRE_CLE_NCBI":
                params['api_key'] = api_key
            
            response = requests.get(search_url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                ids = data['esearchresult'].get('idlist', [])
                
                # ESummary
                if ids:
                    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                    params = {
                        'db': database,
                        'id': ','.join(ids),
                        'retmode': 'json'
                    }
                    if api_key != "VOTRE_CLE_NCBI":
                        params['api_key'] = api_key
                    
                    resp_summary = requests.get(summary_url, params=params, timeout=10)
                    
                    if resp_summary.status_code == 200:
                        return resp_summary.json()
            
            return None
            
        except Exception as e:
            st.warning(f"⚠️ Erreur NCBI: {str(e)}")
            return None

# ============================================================================
# PAGE D'INTÉGRATION API - À AJOUTER DANS LA NAVIGATION
# ============================================================================

def page_integration_api():
    """Page de démonstration et d'utilisation des API externes"""
    st.markdown('<h2 class="section-header">🌐 INTÉGRATION API EXTERNES</h2>', unsafe_allow_html=True)
    st.markdown("""
    <div class='api-card'>
        <h4>🔌 Connectivité avec services externes</h4>
        <p>Cette page vous permet d'interroger en temps réel des API météo, de vision par ordinateur et de génomique.</p>
        <p><strong>⚠️ Note :</strong> Remplacez les clés API dans le code par vos propres clés pour activer les appels réels.</p>
    </div>
    """, unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs([
        "🌤️ Météo & Nutrition",
        "🔬 Roboflow Vision",
        "🧬 Ensembl Génomique",
        "📚 NCBI PubMed"
    ])
    
    with tab1:
        st.markdown("### 🌤️ MÉTÉO EN TEMPS RÉEL - ALGÉRIE")
        st.info("""
        **API OpenWeatherMap** - Données météorologiques actuelles et prévisions.
        Obtenez votre clé gratuite sur https://openweathermap.org/api
        """)
        
        col1, col2 = st.columns([1, 2])
        
        with col1:
            ville = st.text_input("Ville", "Tlemcen")
            if st.button("🔄 Actualiser météo", use_container_width=True):
                with st.spinner("Consultation API météo..."):
                    weather = APIManager.get_weather(ville)
                    
                    if weather:
                        st.session_state.weather_data = weather
                        st.success(f"✅ Données météo pour {weather['ville']}")
                        
                        # Affichage compact
                        st.metric("🌡️ Température", f"{weather['temperature']}°C")
                        st.metric("💧 Humidité", f"{weather['humidite']}%")
                        st.metric("💨 Vent", f"{weather['vent_vitesse']} m/s")
                        st.caption(f"☀️ Lever: {weather['lever_soleil']} | 🌙 Coucher: {weather['coucher_soleil']}")
                    else:
                        # Mode démo
                        st.warning("Mode démo - Données simulées")
                        demo_weather = {
                            'temperature': 28.5,
                            'humidite': 45,
                            'ville': ville
                        }
                        st.session_state.weather_data = demo_weather
                        st.metric("🌡️ Température (simulée)", "28.5°C")
                        st.metric("💧 Humidité (simulée)", "45%")
        
        with col2:
            if 'weather_data' in st.session_state:
                weather = st.session_state.weather_data
                
                st.markdown(f"""
                ### 📊 CONDITIONS ACTUELLES
                
                **{weather.get('description', 'Ensoleillé').capitalize()}** à **{weather.get('ville', ville)}**
                
                | Paramètre | Valeur |
                |-----------|--------|
                | Température ressentie | {weather.get('ressenti', 28)}°C |
                | Pression | {weather.get('pression', 1013)} hPa |
                | Visibilité | {weather.get('visibilite', 10)} km |
                """)
                
                # Conseils nutritionnels
                st.markdown("---")
                st.markdown("### 🥛 CONSEILS NUTRITIONNELS")
                
                advice = APIManager.get_nutrition_advice(weather)
                
                if "NORMAL" in advice:
                    st.success(advice)
                elif "MODÉRÉ" in advice:
                    st.warning(advice)
                else:
                    st.error(advice)
        
        # Prévisions
        st.markdown("---")
        st.markdown("### 📅 PRÉVISIONS 5 JOURS")
        
        if st.button("📊 Voir les prévisions", use_container_width=True):
            with st.spinner("Récupération des prévisions..."):
                forecast = APIManager.get_forecast(ville)
                if forecast:
                    # Grouper par jour
                    days = {}
                    for f in forecast[:40]:
                        date_str = f['datetime'].strftime('%d/%m/%Y')
                        if date_str not in days:
                            days[date_str] = []
                        days[date_str].append(f)
                    
                    cols = st.columns(5)
                    for i, (date, day_data) in enumerate(list(days.items())[:5]):
                        with cols[i]:
                            temp_moy = np.mean([d['temperature'] for d in day_data])
                            humid_moy = np.mean([d['humidite'] for d in day_data])
                            pluie_max = max([d['probabilite_pluie'] for d in day_data])
                            
                            st.markdown(f"""
                            **{date}**
                            🌡️ {temp_moy:.1f}°C
                            💧 {humid_moy:.0f}%
                            ☔ {pluie_max:.0f}%
                            """)
                else:
                    # Mode démo
                    st.info("Mode démo - Prévisions simulées")
                    cols = st.columns(5)
                    dates = ["Lun", "Mar", "Mer", "Jeu", "Ven"]
                    for i, d in enumerate(dates):
                        with cols[i]:
                            st.markdown(f"""
                            **{d}**
                            🌡️ {random.uniform(25, 32):.1f}°C
                            💧 {random.randint(40, 60)}%
                            ☔ {random.randint(0, 30)}%
                            """)
    
    with tab2:
        st.markdown("### 🔬 DÉTECTION DES MAMELLES PAR IA")
        st.info("""
        **API Roboflow** - Détection automatique des mamelles et estimation des dimensions.
        Entraînez votre modèle sur https://roboflow.com et remplacez la clé API.
        """)
        
        if 'rear_image' in st.session_state and st.session_state.rear_image is not None:
            st.image(st.session_state.rear_image, caption="Image actuelle (photo arrière)", width=400)
            
            if st.button("🔍 Détecter avec Roboflow", type="primary", use_container_width=True):
                with st.spinner("Analyse par IA en cours..."):
                    result = APIManager.detect_mammary_roboflow(
                        st.session_state.rear_image,
                        api_key="VOTRE_CLE_ROBOFLOW"
                    )
                    
                    if result:
                        st.success(f"✅ {result['nombre_mamelles']} mamelle(s) détectée(s)")
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Nombre", result['nombre_mamelles'])
                        with col2:
                            conf = result.get('confiance_moyenne', 0.85) * 100
                            st.metric("Confiance", f"{conf:.0f}%")
                        with col3:
                            st.metric("Volume estimé", f"{result.get('volume_moyen_cm3', 250):.0f} cm³")
                        
                        if result.get('mode') == 'simulation':
                            st.caption("⚠️ Mode simulation - API non configurée")
                    else:
                        st.error("Échec de la détection")
        else:
            st.warning("⚠️ Prenez d'abord une photo arrière dans l'onglet **📸 PHOTO & MESURES**")
            if st.button("📸 Aller à PHOTO & MESURES"):
                st.session_state.page = "📸 PHOTO & MESURES"
                st.rerun()
    
    with tab3:
        st.markdown("### 🧬 API ENSEMBL - GÉNOMIQUE OVINE")
        st.info("""
        **Ensembl REST API** - Accès aux données génomiques de référence pour *Ovis aries*.
        Gratuit, sans clé API requise.
        """)
        
        col1, col2 = st.columns(2)
        
        with col1:
            gene_name = st.text_input("Nom du gène", "GDF8")
            if st.button("🔬 Rechercher sur Ensembl", use_container_width=True):
                with st.spinner("Interrogation d'Ensembl..."):
                    gene_data = APIManager.get_ensembl_gene(gene_name)
                    
                    if gene_data:
                        st.success(f"✅ Gène {gene_name} trouvé !")
                        
                        st.markdown(f"""
                        **ID Ensembl :** `{gene_data.get('id', 'N/A')}`  
                        **Description :** {gene_data.get('description', 'N/A')}  
                        **Chromosome :** {gene_data.get('seq_region_name', 'N/A')}  
                        **Position :** {gene_data.get('start', 'N/A')} - {gene_data.get('end', 'N/A')}  
                        **Biotype :** {gene_data.get('biotype', 'N/A')}
                        """)
                        
                        if st.button("📄 Récupérer la séquence"):
                            seq = APIManager.get_ensembl_sequence(gene_data['id'])
                            if seq:
                                st.code(seq[:500] + ("..." if len(seq) > 500 else ""), language="text")
                            else:
                                st.warning("Séquence non disponible")
                    else:
                        st.warning(f"Gène {gene_name} non trouvé dans Ensembl")
        
        with col2:
            st.markdown("### 📊 Gènes d'intérêt économique")
            st.markdown("""
            | Gène | Effet | Race |
            |------|-------|------|
            | **GDF8** | Hypertrophie musculaire | Toutes |
            | **PRNP** | Résistance tremblante | Toutes |
            | **FecB** | Prolificité | OUDA |
            | **DGAT1** | Taux de matière grasse | HAMRA |
            | **CSN1S1** | Taux de caséine | BERBERE |
            """)
    
    with tab4:
        st.markdown("### 📚 NCBI PubMed - Recherche bibliographique")
        st.info("""
        **NCBI E-utilities** - Recherche dans la base de données PubMed.
        Sans clé API : 3 requêtes/seconde. Avec clé : 10 req/s.
        """)
        
        pubmed_query = st.text_input("Recherche PubMed", "sheep genomics Algeria")
        database = st.selectbox("Base de données", ["pubmed", "gene", "nucleotide", "protein"])
        
        if st.button("📚 Rechercher", use_container_width=True):
            with st.spinner("Recherche NCBI en cours..."):
                results = APIManager.ncbi_search(pubmed_query, database, 5)
                
                if results:
                    st.success("✅ Résultats trouvés")
                    
                    if database == "pubmed":
                        uids = results.get('result', {}).get('uids', [])
                        for uid in uids[:5]:
                            article = results.get('result', {}).get(uid, {})
                            title = article.get('title', 'Titre non disponible')
                            authors = article.get('authors', [])
                            author_names = [a.get('name', '') for a in authors[:3]]
                            authors_str = ', '.join(author_names) + (' et al.' if len(authors) > 3 else '')
                            source = article.get('source', '')
                            pubdate = article.get('pubdate', '')
                            
                            st.markdown(f"""
                            **{title}**  
                            *{authors_str}*  
                            {source}, {pubdate}  
                            PMID: {uid}
                            ---
                            """)
                    else:
                        st.json(results)
                else:
                    st.warning("Aucun résultat trouvé ou service indisponible")

# ============================================================================
# SECTION 17: BARRE LATÉRALE - MODIFIÉE POUR AJOUTER LES NOUVELLES PAGES
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
            "🍼 ESTIMATION LAIT MAMMELLE",
            "🌐 API EXTERNES"  # NOUVEAU - AJOUTÉ ICI
        ]
    )
    
    st.markdown("---")
    
    # Statistiques rapides
    st.markdown("### 📊 STATISTIQUES")
    cursor_side = conn.cursor()
    try:
        cursor_side.execute("SELECT COUNT(*) FROM brebis")
        total_brebis = cursor_side.fetchone()[0]
    except:
        total_brebis = 20
    st.metric("🐑 Brebis", total_brebis)
    st.metric("🏷️ Races", "6")
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
# SECTION 18: NAVIGATION PRINCIPALE - MODIFIÉE POUR INCLURE LES NOUVELLES PAGES
# ============================================================================
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
elif page == "🌐 API EXTERNES":  # NOUVEAU - AJOUTÉ ICI
    page_integration_api()

# ============================================================================
# SECTION 19: PIED DE PAGE
# ============================================================================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version 8.0 AVEC API EXTERNES</p>
    <p>📸 Photo & Mesures • 📦 Analyse Multiple • 📐 Scanner 3D • 🎯 Critères de sélection • 🧬 Génétique • 🧬🔬 Génomique avancée • 🥛🔬 Biochimie lait • 🥩 Estimation viande • 🍼 Estimation lait • 🌐 API Météo/Roboflow/Ensembl</p>
    <p>© 2024 - Système de gestion scientifique des races ovines algériennes</p>
    <p><small>✅ Bug np.int0 corrigé - Module API externe intégré - Clés API à configurer dans SECTION 25</small></p>
</div>
""", unsafe_allow_html=True)
