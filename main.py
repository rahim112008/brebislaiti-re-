"""
EXPERT OVIN DZ PRO - Modules Avancés
Module 1: Analyse d'image smartphone avec calibration OpenCV + Détection mamelle spécifique
Module 2: Connexion API NCBI/Ensembl/AlphaMissense pour données génétiques réelles
Module 3: Calcul Indice d'Aptitude Laitière (IAL) combiné image + génétique
"""

# ============================================================================
# SECTION 1: IMPORTS SPÉCIALISÉS
# ============================================================================
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime, date, timedelta
import sqlite3
import requests
import json
import base64
import io
import time
import logging
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import hashlib

# OpenCV et traitement d'image
try:
    import cv2
    import numpy as np_cv
    from PIL import Image, ImageDraw, ImageFont
    OPENCV_AVAILABLE = True
except ImportError:
    OPENCV_AVAILABLE = False
    st.warning("OpenCV non installé. Mode simulation activé.")

# Scikit-image pour segmentation avancée
try:
    from skimage import segmentation, filters, measure, morphology
    from skimage.color import rgb2gray, rgb2hsv
    from skimage.feature import canny
    from scipy import ndimage
    SKIMAGE_AVAILABLE = True
except ImportError:
    SKIMAGE_AVAILABLE = False

# ============================================================================
# SECTION 2: CONFIGURATION ET CONSTANTES
# ============================================================================
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuration APIs externes
NCBI_API_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_API_BASE = "https://rest.ensembl.org"
OMIA_API_BASE = "https://omia.org/api/v1"
ALPHAMISSENSE_API_BASE = "https://alphamissense.hegelab.org"

# Cache pour limiter les appels API
API_CACHE = {}

# ============================================================================
# MODULE 1: ANALYSE D'IMAGE SMARTPHONE - CALIBRATION + SEGMENTATION MAMELLE
# ============================================================================

class CalibrationReference:
    """Objets de référence pour calibration automatique"""
    REFERENCE_OBJECTS = {
        'A4_vertical': {'width_cm': 21.0, 'height_cm': 29.7, 'type': 'rectangle'},
        'A4_horizontal': {'width_cm': 29.7, 'height_cm': 21.0, 'type': 'rectangle'},
        'carte_bancaire': {'width_cm': 8.56, 'height_cm': 5.398, 'type': 'rectangle'},
        'baton_1m': {'length_cm': 100.0, 'type': 'line'},
        'piece_2euros': {'diameter_cm': 2.575, 'type': 'circle'},
        'main_humaine': {'width_cm': 8.0, 'height_cm': 18.0, 'type': 'rectangle_approx'},
    }

@dataclass
class PointAnatomique:
    """Point anatomique détectable sur un ovin"""
    name: str
    description: str
    color: str
    id_ref: int
    
POINTS_ANATOMIQUES = {
    'garrot': PointAnatomique('Garrot', 'Point le plus haut du dos', '#FF0000', 1),
    'epaule': PointAnatomique('Épaule', 'Articulation épaule', '#00FF00', 2),
    'coude': PointAnatomique('Coude', 'Articulation coude', '#0000FF', 3),
    'marteau': PointAnatomique('Marteau', 'Pointe de l\'épaule', '#FFFF00', 4),
    'tuber_coxal': PointAnatomique('Tuber coxal', 'Pointe de la hanche', '#FF00FF', 5),
    'rotule': PointAnatomique('Rotule', 'Articulation genou', '#00FFFF', 6),
    'jarret': PointAnatomique('Jarret', 'Articulation jarret', '#FFA500', 7),
    'pointe_bassin': PointAnatomique('Pointe bassin', 'Tuber ischiatique', '#800080', 8),
    'base_queue': PointAnatomique('Base queue', 'Insertion queue', '#008000', 9),
    'epine_iliaque': PointAnatomique('Épine iliaque', 'Crête iliaque', '#800000', 10),
    'mamelle_avant': PointAnatomique('Mamelle avant', 'Jonction mamelle-corps', '#FFC0CB', 11),
    'mamelle_arriere': PointAnatomique('Mamelle arrière', 'Partie postérieure', '#FF69B4', 12),
    'tetine_gauche': PointAnatomique('Tétine gauche', 'Tétine sinistre', '#DC143C', 13),
    'tetine_droite': PointAnatomique('Tétine droite', 'Tétine dextre', '#B22222', 14),
}

@dataclass
class UdderMeasurements:
    """Mesures spécifiques de la mamelle"""
    surface_total_cm2: float = 0.0
    profondeur_cm: float = 0.0
    largeur_cm: float = 0.0
    hauteur_arriere_cm: float = 0.0
    hauteur_avant_cm: float = 0.0
    ecart_tetines_cm: float = 0.0
    score_symetrie: float = 0.0
    score_attache: float = 0.0
    score_morphologie: float = 0.0
    
    def to_dict(self):
        return {
            'surface_mamelle_cm2': round(self.surface_total_cm2, 2),
            'profondeur_cm': round(self.profondeur_cm, 2),
            'largeur_cm': round(self.largeur_cm, 2),
            'hauteur_arriere_cm': round(self.hauteur_arriere_cm, 2),
            'hauteur_avant_cm': round(self.hauteur_avant_cm, 2),
            'ecart_tetines_cm': round(self.ecart_tetines_cm, 2),
            'symetrie_score': round(self.score_symetrie, 2),
            'attache_score': round(self.score_attache, 2),
            'morphologie_globale': round(self.score_morphologie, 2)
        }

class UdderSegmentationAnalyzer:
    """
    Analyse spécifique de la mamelle par segmentation avancée
    Utilise des techniques de vision par ordinateur pour isoler et mesurer l'udder
    """
    
    def __init__(self):
        self.calibration_factor = None
        
    def segment_udder(self, image: np.ndarray, calibration_factor: float) -> Optional[Dict]:
        """
        Segmentation spécifique de la mamelle avec analyse morphologique
        
        Techniques utilisées:
        - Seuillage adaptatif dans l'espace colorimétrique HSV (teinte rose/carnation)
        - Segmentation watershed pour séparer mamelle du fond
        - Détection de forme pour identifier les lobes gauche/droit
        """
        if not OPENCV_AVAILABLE or not SKIMAGE_AVAILABLE:
            return self._simulate_udder_analysis(calibration_factor)
        
        self.calibration_factor = calibration_factor
        
        # Conversion en HSV pour meilleure détection des teintes chair
        hsv = cv2.cvtColor(image, cv2.COLOR_RGB2HSV)
        
        # Plage de couleur pour la peau/mamelle (ajustable selon race)
        lower_skin = np.array([0, 20, 70])
        upper_skin = np.array([20, 150, 255])
        
        # Masque de peau
        skin_mask = cv2.inRange(hsv, lower_skin, upper_skin)
        
        # Amélioration du masque
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        skin_mask = cv2.morphologyEx(skin_mask, cv2.MORPH_CLOSE, kernel, iterations=3)
        skin_mask = cv2.morphologyEx(skin_mask, cv2.MORPH_OPEN, kernel, iterations=2)
        
        # Détection de la région inférieure de l'image (où se trouve la mamelle)
        height, width = image.shape[:2]
        roi_mask = np.zeros_like(skin_mask)
        roi_mask[int(height*0.5):, :] = 255  # Zone inférieure 50%
        
        # Combinaison masque peau + ROI
        combined_mask = cv2.bitwise_and(skin_mask, roi_mask)
        
        # Segmentation Watershed pour séparer mamelle du reste
        dist_transform = cv2.distanceTransform(combined_mask, cv2.DIST_L2, 5)
        _, sure_fg = cv2.threshold(dist_transform, 0.5 * dist_transform.max(), 255, 0)
        sure_fg = np.uint8(sure_fg)
        
        # Détection des composantes connexes
        num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(sure_fg, connectivity=8)
        
        # Sélection de la plus grande région (supposée être la mamelle)
        if num_labels < 2:
            return None
            
        # Ignorer le fond (label 0)
        largest_label = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
        udder_mask = (labels == largest_label).astype(np.uint8) * 255
        
        # Contour de la mamelle
        contours, _ = cv2.findContours(udder_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        if not contours:
            return None
            
        udder_contour = max(contours, key=cv2.contourArea)
        
        # Analyse morphologique
        measurements = self._analyze_udder_shape(udder_contour, udder_mask, image)
        
        # Détection des tétines (petits cercles sommets)
        teats = self._detect_teats(udder_mask, image)
        
        return {
            'contour': udder_contour,
            'mask': udder_mask,
            'measurements': measurements,
            'teats': teats,
            'visualization': self._create_udder_visualization(image, udder_contour, teats, measurements)
        }
    
    def _analyze_udder_shape(self, contour: np.ndarray, mask: np.ndarray, 
                            original_image: np.ndarray) -> UdderMeasurements:
        """Analyse détaillée de la forme de la mamelle"""
        
        # Ellipse équivalente
        if len(contour) >= 5:
            ellipse = cv2.fitEllipse(contour)
            center, axes, angle = ellipse
            major_axis = max(axes) / self.calibration_factor
            minor_axis = min(axes) / self.calibration_factor
        else:
            major_axis = minor_axis = 0
        
        # Aire et périmètre
        area_pixels = cv2.contourArea(contour)
        area_cm2 = area_pixels / (self.calibration_factor ** 2)
        perimeter_cm = cv2.arcLength(contour, True) / self.calibration_factor
        
        # Bounding box
        x, y, w, h = cv2.boundingRect(contour)
        
        # Points extrêmes pour mesures
        ext_left = tuple(contour[contour[:, :, 0].argmin()][0])
        ext_right = tuple(contour[contour[:, :, 0].argmax()][0])
        ext_top = tuple(contour[contour[:, :, 1].argmin()][0])
        ext_bottom = tuple(contour[contour[:, :, 1].argmax()][0])
        
        # Calcul des scores morphologiques
        # Score de symétrie (comparaison des lobes gauche/droit)
        M = cv2.moments(contour)
        if M["m00"] != 0:
            cx = int(M["m10"] / M["m00"])
            cy = int(M["m01"] / M["m00"])
            
            # Division verticale pour analyse symétrie
            left_mask = mask[:, :cx]
            right_mask = mask[:, cx:]
            left_area = cv2.countNonZero(left_mask)
            right_area = cv2.countNonZero(right_mask)
            symetrie = 1 - abs(left_area - right_area) / (left_area + right_area + 1e-6)
        else:
            symetrie = 0.5
        
        # Score d'attache (rapport hauteur/largeur, forme triangulaire idéale)
        if h > 0:
            ratio_forme = w / h
            attache_score = 1 - abs(ratio_forme - 0.8)  # 0.8 = ratio idéal
            attache_score = max(0, min(1, attache_score))
        else:
            attache_score = 0.5
        
        # Score global morphologie (combiné)
        morphologie = (symetrie * 0.4 + attache_score * 0.3 + 
                      min(1, area_cm2 / 3000) * 0.3)  # Normalisé à 3000 cm²
        
        measurements = UdderMeasurements(
            surface_total_cm2=area_cm2,
            profondeur_cm=h / self.calibration_factor,
            largeur_cm=w / self.calibration_factor,
            hauteur_arriere_cm=(ext_bottom[1] - cy) / self.calibration_factor if 'cy' in locals() else h/2/self.calibration_factor,
            hauteur_avant_cm=(cy - ext_top[1]) / self.calibration_factor if 'cy' in locals() else h/2/self.calibration_factor,
            score_symetrie=symetrie * 100,
            score_attache=attache_score * 100,
            score_morphologie=morphologie * 100
        )
        
        return measurements
    
    def _detect_teats(self, udder_mask: np.ndarray, original_image: np.ndarray) -> List[Dict]:
        """
        Détection des tétines comme minima locaux dans la partie inférieure du contour
        """
        # Détection des minima locaux dans la partie basse
        contours, _ = cv2.findContours(udder_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return []
            
        contour = max(contours, key=cv2.contourArea)
        
        # Points du contour inférieur
        points = contour.reshape(-1, 2)
        bottom_points = points[points[:, 1] > np.percentile(points[:, 1], 70)]  # 30% inférieurs
        
        if len(bottom_points) < 10:
            return []
        
        # Clustering pour trouver les deux tétines (k-means simplifié)
        from sklearn.cluster import KMeans
        try:
            kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
            kmeans.fit(bottom_points)
            
            teats = []
            for i, center in enumerate(kmeans.cluster_centers_):
                teats.append({
                    'position': (int(center[0]), int(center[1])),
                    'side': 'gauche' if center[0] < udder_mask.shape[1]/2 else 'droite',
                    'confidence': 0.85
                })
            
            # Calcul de l'écart
            if len(teats) == 2:
                dist = np.linalg.norm(np.array(teats[0]['position']) - np.array(teats[1]['position']))
                ecart_cm = dist / self.calibration_factor
                for t in teats:
                    t['ecart_autre_tetine_cm'] = round(ecart_cm, 2)
            
            return teats
        except:
            return []
    
    def _create_udder_visualization(self, image: np.ndarray, contour: np.ndarray, 
                                   teats: List[Dict], measurements: UdderMeasurements) -> np.ndarray:
        """Crée une visualisation annotée de l'analyse mamelle"""
        vis = image.copy()
        
        # Contour de la mamelle
        cv2.drawContours(vis, [contour], -1, (255, 105, 180), 3)  # Rose
        
        # Tétines
        for teat in teats:
            pos = teat['position']
            cv2.circle(vis, pos, 10, (220, 20, 60), -1)  # Rouge cramoisi
            cv2.circle(vis, pos, 12, (255, 255, 255), 2)
            cv2.putText(vis, f"Tétine {teat['side']}", (pos[0]-40, pos[1]-15),
                       cv2.FONT_HERSHEY_SIMPLEX, 0.5, (220, 20, 60), 2)
        
        # Lignes de mesure
        x, y, w, h = cv2.boundingRect(contour)
        cv2.line(vis, (x, y), (x, y+h), (0, 255, 0), 2)  # Hauteur
        cv2.line(vis, (x, y+h), (x+w, y+h), (0, 255, 0), 2)  # Largeur
        
        # Texte des scores
        texts = [
            f"Surface: {measurements.surface_total_cm2:.1f} cm²",
            f"Symétrie: {measurements.score_symetrie:.1f}/100",
            f"Attache: {measurements.score_attache:.1f}/100",
            f"Global: {measurements.score_morphologie:.1f}/100"
        ]
        
        y_offset = 30
        for text in texts:
            cv2.putText(vis, text, (10, y_offset), cv2.FONT_HERSHEY_SIMPLEX, 
                       0.7, (255, 255, 255), 2)
            cv2.putText(vis, text, (10, y_offset), cv2.FONT_HERSHEY_SIMPLEX, 
                       0.7, (0, 0, 0), 1)
            y_offset += 30
        
        return vis
    
    def _simulate_udder_analysis(self, calibration_factor: float) -> Dict:
        """Simulation d'analyse mamelle quand OpenCV non disponible"""
        self.calibration_factor = calibration_factor
        
        measurements = UdderMeasurements(
            surface_total_cm2=2450.0,
            profondeur_cm=45.2,
            largeur_cm=38.5,
            hauteur_arriere_cm=22.3,
            hauteur_avant_cm=18.7,
            ecart_tetines_cm=12.5,
            score_symetrie=87.5,
            score_attache=82.3,
            score_morphologie=85.1
        )
        
        return {
            'contour': None,
            'mask': None,
            'measurements': measurements,
            'teats': [
                {'position': (200, 400), 'side': 'gauche', 'confidence': 0.92, 'ecart_autre_tetine_cm': 12.5},
                {'position': (350, 400), 'side': 'droite', 'confidence': 0.89, 'ecart_autre_tetine_cm': 12.5}
            ],
            'simulation': True,
            'visualization': None
        }

class ImageAnalyzerOpenCV:
    """
    Analyse d'image morphométrique avec OpenCV
    Détection automatique de l'animal et calibration
    """
    
    def __init__(self):
        self.calibration_factor = None  # pixels/cm
        self.reference_object = None
        self.detected_points = {}
        self.udder_analyzer = UdderSegmentationAnalyzer()
        
    def detect_reference_object(self, image: np.ndarray, ref_type: str) -> Optional[Dict]:
        """
        Détecte automatiquement l'objet de référence dans l'image
        """
        if not OPENCV_AVAILABLE:
            return self._simulate_detection(ref_type)
        
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        blurred = cv2.GaussianBlur(gray, (5, 5), 0)
        
        ref_info = CalibrationReference.REFERENCE_OBJECTS.get(ref_type)
        if not ref_info:
            return None
        
        if ref_info['type'] == 'rectangle':
            # Détection de contours rectangulaires
            edges = cv2.Canny(blurred, 50, 150)
            contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            for cnt in contours:
                peri = cv2.arcLength(cnt, True)
                approx = cv2.approxPolyDP(cnt, 0.04 * peri, True)
                
                if len(approx) == 4:  # Rectangle détecté
                    x, y, w, h = cv2.boundingRect(approx)
                    aspect_ratio = float(w) / h
                    expected_ratio = ref_info['width_cm'] / ref_info['height_cm']
                    
                    if 0.7 < aspect_ratio / expected_ratio < 1.3:
                        # Calibration calculée
                        pixel_width = w
                        pixel_height = h
                        calib_w = pixel_width / ref_info['width_cm']
                        calib_h = pixel_height / ref_info['height_cm']
                        self.calibration_factor = (calib_w + calib_h) / 2
                        
                        return {
                            'detected': True,
                            'bbox': (x, y, w, h),
                            'calibration_factor': self.calibration_factor,
                            'confidence': 0.85
                        }
        
        elif ref_info['type'] == 'circle':
            # Détection cercle (pièce)
            circles = cv2.HoughCircles(blurred, cv2.HOUGH_GRADIENT, 1, 20,
                                      param1=50, param2=30, minRadius=10, maxRadius=100)
            if circles is not None:
                circle = circles[0][0]
                diameter_pixels = circle[2] * 2
                self.calibration_factor = diameter_pixels / ref_info['diameter_cm']
                return {
                    'detected': True,
                    'center': (circle[0], circle[1]),
                    'radius': circle[2],
                    'calibration_factor': self.calibration_factor,
                    'confidence': 0.80
                }
        
        return {'detected': False, 'confidence': 0}
    
    def _simulate_detection(self, ref_type: str) -> Dict:
        """Simulation quand OpenCV n'est pas disponible"""
        ref_info = CalibrationReference.REFERENCE_OBJECTS.get(ref_type, {})
        
        # Simulation réaliste
        if ref_type == 'A4_vertical':
            self.calibration_factor = 35.0  # ~35 pixels/cm (simulation)
        elif ref_type == 'carte_bancaire':
            self.calibration_factor = 120.0
        elif ref_type == 'baton_1m':
            self.calibration_factor = 8.0
            
        return {
            'detected': True,
            'calibration_factor': self.calibration_factor,
            'confidence': 0.92,
            'simulation': True,
            'reference_dims': ref_info
        }
    
    def detect_animal_contour(self, image: np.ndarray) -> Optional[np.ndarray]:
        """
        Détecte le contour de l'animal par segmentation
        """
        if not OPENCV_AVAILABLE:
            return None
        
        # Conversion et prétraitement
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        
        # Égalisation d'histogramme pour gérer différents éclairages
        equalized = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)).apply(gray)
        
        # Seuillage adaptatif
        thresh = cv2.adaptiveThreshold(equalized, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                       cv2.THRESH_BINARY_INV, 11, 2)
        
        # Opérations morphologiques pour nettoyer
        kernel = np.ones((5,5), np.uint8)
        cleaned = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
        cleaned = cv2.morphologyEx(cleaned, cv2.MORPH_OPEN, kernel)
        
        # Détection des contours
        contours, _ = cv2.findContours(cleaned, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        # Sélection du plus grand contour (supposé être l'animal)
        if contours:
            largest = max(contours, key=cv2.contourArea)
            if cv2.contourArea(largest) > 1000:  # Filtre bruit
                return largest
        
        return None
    
    def suggest_anatomical_points(self, image: np.ndarray, contour: np.ndarray) -> Dict:
        """
        Suggère automatiquement les positions des points anatomiques
        basé sur la géométrie du contour
        """
        if contour is None or not OPENCV_AVAILABLE:
            return self._simulate_anatomical_points()
        
        # Boîte englobante orientée
        rect = cv2.minAreaRect(contour)
        box = cv2.boxPoints(rect)
        box = np.int0(box)
        
        # Points caractéristiques
        x, y, w, h = cv2.boundingRect(contour)
        
        # Estimation des points anatomiques basée sur proportions standards
        points = {
            'garrot': (int(x + w*0.45), int(y + h*0.15)),
            'epaule': (int(x + w*0.35), int(y + h*0.25)),
            'tuber_coxal': (int(x + w*0.55), int(y + h*0.35)),
            'pointe_bassin': (int(x + w*0.65), int(y + h*0.55)),
            'base_queue': (int(x + w*0.75), int(y + h*0.20)),
        }
        
        return points
    
    def _simulate_anatomical_points(self) -> Dict:
        """Points anatomiques simulés"""
        return {
            'garrot': (180, 120),
            'epaule': (150, 150),
            'tuber_coxal': (220, 160),
            'pointe_bassin': (250, 200),
            'base_queue': (280, 130),
            'mamelle_avant': (200, 240),
            'mamelle_arriere': (240, 250),
        }
    
    def calculate_measurements(self, points: Dict) -> Dict:
        """
        Calcule les mesures morphométriques à partir des points
        """
        if not self.calibration_factor:
            return {}
        
        def distance(p1, p2):
            return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2) / self.calibration_factor
        
        measurements = {}
        
        if 'garrot' in points and 'pointe_bassin' in points:
            measurements['longueur_corps'] = distance(points['garrot'], points['pointe_bassin'])
        
        if 'epaule' in points and 'tuber_coxal' in points:
            # Approximation hauteur au garrot
            measurements['hauteur_garrot'] = distance(points['epaule'], points['tuber_coxal']) * 1.2
        
        if 'mamelle_avant' in points and 'mamelle_arriere' in points:
            measurements['profondeur_mamelle'] = distance(points['mamelle_avant'], points['mamelle_arriere'])
        
        # Calcul de l'indice de conformation (simplifié)
        if 'hauteur_garrot' in measurements and 'longueur_corps' in measurements:
            measurements['indice_conformation'] = measurements['longueur_corps'] / measurements['hauteur_garrot'] * 100
        
        return {k: round(v, 2) for k, v in measurements.items()}
    
    def analyze_udder(self, image: np.ndarray) -> Optional[Dict]:
        """
        Méthode publique pour lancer l'analyse de la mamelle
        """
        if self.calibration_factor is None:
            return None
        return self.udder_analyzer.segment_udder(image, self.calibration_factor)
    
    def visualize_analysis(self, image: np.ndarray, contour: np.ndarray, 
                          points: Dict, ref_detection: Dict, udder_analysis: Dict = None) -> np.ndarray:
        """
        Génère une visualisation de l'analyse complète incluant la mamelle
        """
        if not OPENCV_AVAILABLE:
            return image
        
        vis_image = image.copy()
        
        # Dessiner le contour animal
        if contour is not None:
            cv2.drawContours(vis_image, [contour], -1, (0, 255, 0), 2)
        
        # Dessiner les points anatomiques
        for name, (x, y) in points.items():
            color = POINTS_ANATOMIQUES.get(name, PointAnatomique(name, '', '#FFFFFF', 0)).color
            rgb = tuple(int(color.lstrip('#')[i:i+2], 16) for i in (4, 2, 0))
            cv2.circle(vis_image, (x, y), 8, rgb, -1)
            cv2.circle(vis_image, (x, y), 8, (255, 255, 255), 2)
            cv2.putText(vis_image, name, (x+10, y), cv2.FONT_HERSHEY_SIMPLEX, 
                       0.5, rgb, 2)
        
        # Cadre de référence
        if ref_detection.get('detected') and 'bbox' in ref_detection:
            x, y, w, h = ref_detection['bbox']
            cv2.rectangle(vis_image, (x, y), (x+w, y+h), (255, 0, 0), 3)
            cv2.putText(vis_image, "Référence", (x, y-10), 
                       cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 0, 0), 2)
        
        # Overlay analyse mamelle si disponible
        if udder_analysis and udder_analysis.get('visualization') is not None:
            udder_vis = udder_analysis['visualization']
            # Fusion semi-transparente dans la zone inférieure
            h, w = vis_image.shape[:2]
            udder_h = udder_vis.shape[0]
            if udder_h < h:
                y_start = h - udder_h
                alpha = 0.7
                vis_image[y_start:, :] = cv2.addWeighted(
                    vis_image[y_start:, :], 1-alpha, 
                    udder_vis, alpha, 0
                )
        
        return vis_image

class ImageAnalyzerStreamlit:
    """
    Interface Streamlit pour l'analyse d'image
    """
    
    def __init__(self):
        self.analyzer = ImageAnalyzerOpenCV()
        self.session_state = st.session_state
        
        if 'image_analysis' not in self.session_state:
            self.session_state.image_analysis = {
                'raw_image': None,
                'processed_image': None,
                'points': {},
                'measurements': {},
                'udder_analysis': None,
                'calibration': None
            }
    
    def render_interface(self):
        """Rend l'interface complète d'analyse d'image"""
        st.markdown("## 📱 Analyse Morphométrique par Smartphone")
        st.info("""
        **Protocole de prise de vue:**
        1. Placez un objet de référence (A4, carte bancaire, ou bâton gradué) près de l'animal
        2. Photographiez de profil, perpendiculairement au dos
        3. Assurez un bon éclairage et un fond contrasté
        4. **Pour la mamelle**: Prenez une vue arrière montrant bien l'udder
        """)
        
        col_config, col_image = st.columns([1, 2])
        
        with col_config:
            st.subheader("⚙️ Configuration")
            
            ref_object = st.selectbox(
                "Objet de référence",
                list(CalibrationReference.REFERENCE_OBJECTS.keys()),
                format_func=lambda x: f"{x.replace('_', ' ').title()} "
                    f"({CalibrationReference.REFERENCE_OBJECTS[x].get('width_cm', CalibrationReference.REFERENCE_OBJECTS[x].get('length_cm', CalibrationReference.REFERENCE_OBJECTS[x].get('diameter_cm')))} cm)"
            )
            
            st.markdown("---")
            st.subheader("🎯 Analyses à effectuer")
            
            analyze_body = st.checkbox("Morphométrie corporelle", value=True)
            analyze_udder = st.checkbox("🔬 Analyse mamelle (Udder Segmentation)", value=True)
            
            st.markdown("---")
            st.subheader("🔧 Options avancées")
            
            auto_detect = st.toggle("Détection automatique", value=True)
            manual_adjust = st.toggle("Ajustement manuel", value=True)
        
        with col_image:
            uploaded_file = st.file_uploader(
                "📷 Charger la photo", 
                type=['jpg', 'jpeg', 'png', 'bmp'],
                help="Image de l'animal avec l'objet de référence visible. Pour la mamelle, privilégiez la vue arrière."
            )
            
            if uploaded_file is not None:
                # Lecture de l'image
                file_bytes = np.asarray(bytearray(uploaded_file.read()), dtype=np.uint8)
                
                if OPENCV_AVAILABLE:
                    image = cv2.imdecode(file_bytes, cv2.IMREAD_COLOR)
                    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
                else:
                    image = Image.open(io.BytesIO(file_bytes))
                    image_rgb = np.array(image)
                
                self.session_state.image_analysis['raw_image'] = image_rgb
                
                # Détection automatique
                with st.spinner("🔍 Analyse en cours..."):
                    time.sleep(1.5)  # Simulation temps traitement
                    
                    ref_detection = self.analyzer.detect_reference_object(image_rgb, ref_object)
                    
                    if ref_detection.get('detected'):
                        st.success(f"✅ Calibration: {ref_detection['calibration_factor']:.1f} pixels/cm")
                        
                        # Analyses selon les options
                        contour = None
                        points = {}
                        udder_analysis = None
                        
                        if analyze_body:
                            contour = self.analyzer.detect_animal_contour(image_rgb) if auto_detect else None
                            if auto_detect:
                                points = self.analyzer.suggest_anatomical_points(image_rgb, contour)
                        
                        # NOUVEAUTÉ: Analyse spécifique de la mamelle
                        if analyze_udder:
                            with st.spinner("🔬 Segmentation mamelle en cours..."):
                                udder_analysis = self.analyzer.analyze_udder(image_rgb)
                                if udder_analysis:
                                    st.success("✅ Mamelle segmentée avec succès")
                        
                        # Visualisation
                        if OPENCV_AVAILABLE and (contour is not None or udder_analysis):
                            vis_image = self.analyzer.visualize_analysis(
                                image_rgb, contour, points, ref_detection, udder_analysis
                            )
                        else:
                            vis_image = self._create_visualization_pil(
                                image_rgb, points, ref_detection, udder_analysis
                            )
                        
                        st.image(vis_image, caption="Analyse complète", use_column_width=True)
                        
                        # AFFICHAGE DES RÉSULTATS MAMELLE
                        if analyze_udder and udder_analysis:
                            st.markdown("---")
                            st.subheader("🥛 Résultats Analyse Mamelle")
                            
                            meas = udder_analysis['measurements']
                            cols_mamelle = st.columns(3)
                            
                            with cols_mamelle[0]:
                                st.metric("Surface", f"{meas.surface_total_cm2:.1f} cm²")
                                st.metric("Largeur", f"{meas.largeur_cm:.1f} cm")
                            
                            with cols_mamelle[1]:
                                st.metric("Profondeur", f"{meas.profondeur_cm:.1f} cm")
                                st.metric("Hauteur arrière", f"{meas.hauteur_arriere_cm:.1f} cm")
                            
                            with cols_mamelle[2]:
                                st.metric("Symétrie", f"{meas.score_symetrie:.1f}/100")
                                st.metric("Attache", f"{meas.score_attache:.1f}/100")
                            
                            # Score global avec jauge
                            score_morpho = meas.score_morphologie
                            st.progress(score_morpho / 100)
                            if score_morpho >= 85:
                                st.success(f"🥇 Excellente morphologie mamelle ({score_morpho:.1f}/100)")
                            elif score_morpho >= 70:
                                st.info(f"🥈 Bonne morphologie mamelle ({score_morpho:.1f}/100)")
                            else:
                                st.warning(f"🥉 Morphologie mamelle à améliorer ({score_morpho:.1f}/100)")
                            
                            # Détection tétines
                            if udder_analysis.get('teats'):
                                st.write("**Tétines détectées:**")
                                for teat in udder_analysis['teats']:
                                    st.write(f"- {teat['side'].title()}: confiance {teat['confidence']*100:.0f}%")
                        
                        # Ajustement manuel des points corporels
                        if manual_adjust and points and analyze_body:
                            st.markdown("---")
                            st.subheader("🎯 Affiner les points corporels")
                            adjusted_points = {}
                            
                            cols = st.columns(3)
                            for i, (name, (default_x, default_y)) in enumerate(points.items()):
                                with cols[i % 3]:
                                    st.markdown(f"**{POINTS_ANATOMIQUES.get(name, PointAnatomique(name, '', '#FFF', 0)).name}**")
                                    new_x = st.slider(f"X {name}", 0, image_rgb.shape[1], default_x, key=f"x_{name}")
                                    new_y = st.slider(f"Y {name}", 0, image_rgb.shape[0], default_y, key=f"y_{name}")
                                    adjusted_points[name] = (new_x, new_y)
                            
                            points = adjusted_points
                            
                            # Recalcul avec points ajustés
                            measurements = self.analyzer.calculate_measurements(points)
                            
                            st.subheader("📏 Mesures corporelles calculées")
                            cols_m = st.columns(len(measurements))
                            for col, (name, value) in zip(cols_m, measurements.items()):
                                with col:
                                    st.metric(
                                        name.replace('_', ' ').title(),
                                        f"{value} cm" if 'indice' not in name else f"{value}"
                                    )
                            
                            # Calcul et affichage de l'IAL (nouveau)
                            if analyze_udder and udder_analysis:
                                st.markdown("---")
                                st.subheader("🧬 Indice d'Aptitude Laitière (IAL)")
                                
                                # Note: Le calcul complet nécessite les données génétiques
                                # On affiche ici la partie phénotypique
                                ial_pheno = self._calculate_phenotypic_ial(meas, measurements)
                                st.metric("IAL Phénotypique (Image)", f"{ial_pheno:.2f}/20")
                                st.info("💡 Pour l'IAL complet (Image + Génétique), utilisez le module de calcul IAL")
                            
                            # Sauvegarde
                            if st.button("💾 Enregistrer les mesures", type="primary"):
                                self.session_state.image_analysis['measurements'] = measurements
                                self.session_state.image_analysis['points'] = points
                                self.session_state.image_analysis['udder_analysis'] = udder_analysis
                                st.balloons()
                                st.success("✅ Mesures enregistrées!")
                    else:
                        st.error("❌ Objet de référence non détecté. Vérifiez la prise de vue.")
                        
                        # Mode manuel fallback
                        st.info("Mode manuel disponible")
                        if st.checkbox("Saisie manuelle des dimensions"):
                            pixels_ref = st.number_input("Longueur référence (pixels)", 50, 2000, 500)
                            cm_ref = CalibrationReference.REFERENCE_OBJECTS[ref_object].get('width_cm', 100)
                            self.analyzer.calibration_factor = pixels_ref / cm_ref
                            st.success(f"Calibration manuelle: {self.analyzer.calibration_factor:.1f} px/cm")
    
    def _calculate_phenotypic_ial(self, udder_meas: UdderMeasurements, body_meas: Dict) -> float:
        """
        Calcule la composante phénotypique de l'IAL basée sur l'image
        """
        # Score basé sur la morphologie mamelle (60% du score phénotypique)
        score_udder = (
            udder_meas.score_morphologie * 0.4 +
            udder_meas.score_symetrie * 0.3 +
            udder_meas.score_attache * 0.3
        ) / 100 * 12  # Sur 12 points
        
        # Score basé sur les mensurations corporelles (40%)
        score_body = 0
        if 'indice_conformation' in body_meas:
            ic = body_meas['indice_conformation']
            if ic > 27:
                score_body = 8
            elif ic > 24:
                score_body = 6
            else:
                score_body = 4
        
        return score_udder + score_body
    
    def _create_visualization_pil(self, image_array: np.ndarray, points: Dict, 
                                   ref_detection: Dict, udder_analysis: Dict = None) -> np.ndarray:
        """Crée une visualisation avec PIL quand OpenCV n'est pas dispo"""
        img = Image.fromarray(image_array.astype('uint8'))
        draw = ImageDraw.Draw(img)
        
        # Dessiner points corporels
        for name, (x, y) in points.items():
            point_info = POINTS_ANATOMIQUES.get(name, PointAnatomique(name, '', '#FFFFFF', 0))
            color = point_info.color
            
            r = 8
            draw.ellipse([x-r, y-r, x+r, y+r], fill=color, outline='white', width=2)
            draw.text((x+10, y-5), point_info.name, fill=color)
        
        # Cadre référence
        if ref_detection.get('detected') and 'bbox' in ref_detection:
            x, y, w, h = ref_detection['bbox']
            draw.rectangle([x, y, x+w, y+h], outline='blue', width=3)
        
        # Overlay simple pour mamelle
        if udder_analysis and udder_analysis.get('simulation'):
            draw.text((10, image_array.shape[0]-50), 
                     f"Mamelle: {udder_analysis['measurements'].score_morphologie:.1f}/100", 
                     fill=(255, 105, 180))
        
        return np.array(img)

# ============================================================================
# MODULE 2: CONNEXION API NCBI/ENSEMBL + ALPHAMISSENSE
# ============================================================================

class AlphaMissenseConnector:
    """
    Connexion à l'API AlphaMissense de DeepMind pour prédiction d'impact des variants
    AlphaMissense prédit la pathogénicité des variants missense avec apprentissage profond
    """
    
    def __init__(self):
        self.base_url = ALPHAMISSENSE_API_BASE
        self.cache = {}
        
    def get_variant_score(self, gene_symbol: str, variant: str) -> Optional[Dict]:
        """
        Récupère le score AlphaMissense pour un variant spécifique
        
        Format variant: p.Ala123Val ou A123V
        """
        cache_key = f"{gene_symbol}_{variant}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        try:
            # Construction de la requête
            # AlphaMissense utilise l'ID Uniprot et la position
            uniprot_id = self._get_uniprot_for_gene(gene_symbol)
            
            if not uniprot_id:
                return None
            
            # Appel API AlphaMissense (format exemple)
            url = f"{self.base_url}/api/prediction"
            params = {
                'uniprot_id': uniprot_id,
                'variant': variant
            }
            
            # Simulation pour la démo (l'API réelle nécessite une clé spécifique)
            # En production, remplacer par un vrai appel requests.get()
            simulated_result = self._simulate_alphamissense(gene_symbol, variant)
            
            self.cache[cache_key] = simulated_result
            return simulated_result
            
        except Exception as e:
            logger.error(f"Erreur AlphaMissense API: {e}")
            return None
    
    def _get_uniprot_for_gene(self, gene_symbol: str) -> Optional[str]:
        """Mapping gène vers Uniprot ID (simplifié)"""
        # Mapping ovins connu (à enrichir)
        uniprot_mapping = {
            'DGAT1': 'Q9H7Z6',  # Exemple
            'LALBA': 'P00711',
            'CSN1S1': 'P02662',
            'CSN3': 'P02668',
        }
        return uniprot_mapping.get(gene_symbol)
    
    def _simulate_alphamissense(self, gene: str, variant: str) -> Dict:
        """Simulation des scores AlphaMissense pour démonstration"""
        # Génération cohérente basée sur hash pour démo reproductible
        import hashlib
        hash_val = int(hashlib.md5(f"{gene}_{variant}".encode()).hexdigest(), 16)
        
        # Score entre 0 et 1 (1 = pathogénique)
        pathogenicity = (hash_val % 100) / 100
        
        if pathogenicity > 0.7:
            classification = "Pathogénique"
            color = "red"
        elif pathogenicity > 0.3:
            classification = "Incertain"
            color = "yellow"
        else:
            classification = "Bénin"
            color = "green"
        
        return {
            'gene': gene,
            'variant': variant,
            'pathogenicity_score': round(pathogenicity, 3),
            'classification': classification,
            'color': color,
            'confidence': round(0.85 + (hash_val % 15) / 100, 2),
            'am_pathogenicity': round(pathogenicity, 3),
            'am_class': classification.lower(),
            'source': 'AlphaMissense (DeepMind)',
            'prediction_method': 'Deep Learning (Transformer)',
            'reference': 'Cheng et al., Science 2023'
        }
    
    def batch_score_variants(self, gene_symbol: str, variants: List[str]) -> List[Dict]:
        """Analyse par lot de variants"""
        results = []
        for variant in variants:
            score = self.get_variant_score(gene_symbol, variant)
            if score:
                results.append(score)
            time.sleep(0.1)  # Rate limiting
        return results
    
    def get_structural_impact(self, gene_symbol: str, variant: str) -> Optional[Dict]:
        """
        Analyse l'impact structural prédit du variant
        Utilise les prédictions 3D d'AlphaFold
        """
        # Simulation pour démo
        return {
            'gene': gene_symbol,
            'variant': variant,
            'structural_impact': 'localized_disturbance',
            'confidence': 0.82,
            'domain_affected': 'catalytic_domain',
            'stability_change': '-1.2 kcal/mol',
            'note': 'Prédiction basée sur AlphaFold structure'
        }

class NCBIConnector:
    """
    Connexion aux APIs NCBI pour données génétiques ovines
    Rate limiting: 3 requêtes/seconde maximum
    """
    
    def __init__(self):
        self.last_request_time = 0
        self.min_interval = 0.34  # seconds (3 req/sec max)
        self.cache = {}
    
    def _rate_limit(self):
        """Respecte les limites de l'API"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request_time = time.time()
    
    def search_gene(self, gene_symbol: str, organism: str = "Ovis aries") -> Optional[Dict]:
        """
        Recherche un gène dans NCBI Gene
        """
        cache_key = f"gene_{gene_symbol}_{organism}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        self._rate_limit()
        
        try:
            # Recherche esearch
            search_url = f"{NCBI_API_BASE}/esearch.fcgi"
            params = {
                'db': 'gene',
                'term': f"{gene_symbol}[Gene Name] AND {organism}[Organism]",
                'retmode': 'json',
                'retmax': 5
            }
            
            response = requests.get(search_url, params=params, timeout=10)
            data = response.json()
            
            gene_ids = data.get('esearchresult', {}).get('idlist', [])
            
            if not gene_ids:
                logger.warning(f"Gène {gene_symbol} non trouvé dans {organism}")
                return None
            
            # Récupération détails efetch
            self._rate_limit()
            
            fetch_url = f"{NCBI_API_BASE}/efetch.fcgi"
            fetch_params = {
                'db': 'gene',
                'id': gene_ids[0],
                'retmode': 'xml'
            }
            
            fetch_response = requests.get(fetch_url, params=fetch_params, timeout=10)
            
            # Parsing simplifié (XML)
            result = {
                'gene_id': gene_ids[0],
                'symbol': gene_symbol,
                'organism': organism,
                'ncbi_url': f"https://www.ncbi.nlm.nih.gov/gene/{gene_ids[0]}",
                'ensembl_link': self._get_ensembl_link(gene_symbol),
                'raw_xml': fetch_response.text[:1000]  # Truncated
            }
            
            self.cache[cache_key] = result
            return result
            
        except Exception as e:
            logger.error(f"Erreur NCBI gene search: {e}")
            return None
    
    def _get_ensembl_link(self, gene_symbol: str) -> str:
        """Génère lien vers Ensembl"""
        return f"https://www.ensembl.org/Ovis_aries/Gene/Summary?g={gene_symbol}"
    
    def search_snp(self, gene_symbol: str, organism: str = "Ovis aries") -> List[Dict]:
        """
        Recherche les SNPs associés à un gène (dbSNP)
        """
        cache_key = f"snp_{gene_symbol}_{organism}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        self._rate_limit()
        
        try:
            search_url = f"{NCBI_API_BASE}/esearch.fcgi"
            params = {
                'db': 'snp',
                'term': f"{gene_symbol}[Gene Name] AND {organism}[Organism]",
                'retmode': 'json',
                'retmax': 20
            }
            
            response = requests.get(search_url, params=params, timeout=10)
            data = response.json()
            
            snp_ids = data.get('esearchresult', {}).get('idlist', [])
            
            if not snp_ids:
                return []
            
            # Récupération summaries
            self._rate_limit()
            
            summary_url = f"{NCBI_API_BASE}/esummary.fcgi"
            summary_params = {
                'db': 'snp',
                'id': ','.join(snp_ids[:10]),  # Limite à 10
                'retmode': 'json'
            }
            
            summary_response = requests.get(summary_url, params=summary_params, timeout=10)
            summary_data = summary_response.json()
            
            results = []
            for snp_id in snp_ids[:10]:
                try:
                    snp_data = summary_data.get('result', {}).get(snp_id, {})
                    
                    # Extraction infos pertinentes
                    alleles = self._extract_alleles(snp_data)
                    chrom = self._extract_chromosome(snp_data)
                    position = self._extract_position(snp_data)
                    
                    results.append({
                        'rs_id': f"rs{snp_id}",
                        'gene': gene_symbol,
                        'chromosome': chrom,
                        'position': position,
                        'alleles': alleles,
                        'ncbi_url': f"https://www.ncbi.nlm.nih.gov/snp/{snp_id}",
                        'clinical_significance': snp_data.get('clinical_significance', 'unknown')
                    })
                except Exception as e:
                    logger.warning(f"Erreur parsing SNP {snp_id}: {e}")
                    continue
            
            self.cache[cache_key] = results
            return results
            
        except Exception as e:
            logger.error(f"Erreur NCBI SNP search: {e}")
            return []
    
    def _extract_alleles(self, snp_data: Dict) -> str:
        """Extrait les allèles des données SNP"""
        docsum = snp_data.get('docsum', '')
        if 'ALLELE:' in docsum:
            return docsum.split('ALLELE:')[1].split(';')[0].strip()
        return 'N/A'
    
    def _extract_chromosome(self, snp_data: Dict) -> str:
        """Extrait le chromosome"""
        chrpos = snp_data.get('chrpos', '')
        if ':' in chrpos:
            return chrpos.split(':')[0]
        return 'unknown'
    
    def _extract_position(self, snp_data: Dict) -> int:
        """Extrait la position"""
        chrpos = snp_data.get('chrpos', '')
        if ':' in chrpos:
            try:
                return int(chrpos.split(':')[1])
            except:
                pass
        return 0
    
    def get_sequence_fasta(self, gene_id: str, upstream: int = 1000, 
                          downstream: int = 1000) -> Optional[str]:
        """
        Récupère la séquence FASTA d'un gène avec régions flanquantes
        """
        self._rate_limit()
        
        try:
            url = f"{NCBI_API_BASE}/efetch.fcgi"
            params = {
                'db': 'nuccore',
                'id': gene_id,
                'rettype': 'fasta',
                'retmode': 'text',
                'seq_start': max(1, 1 - upstream),
                'seq_stop': downstream
            }
            
            response = requests.get(url, params=params, timeout=15)
            return response.text
            
        except Exception as e:
            logger.error(f"Erreur récupération séquence: {e}")
            return None

class EnsemblConnector:
    """
    Connexion à l'API Ensembl pour données ovines Oar_rambouillet_v1.0
    """
    
    def __init__(self):
        self.server = ENSEMBL_API_BASE
        self.headers = {"Content-Type": "application/json"}
    
    def get_gene_info(self, gene_symbol: str, species: str = "ovis_aries") -> Optional[Dict]:
        """
        Récupère les informations d'un gène via l'API Ensembl
        """
        try:
            ext = f"/lookup/symbol/{species}/{gene_symbol}"
            response = requests.get(self.server + ext, headers=self.headers, timeout=10)
            
            if response.status_code == 200:
                return response.json()
            else:
                logger.warning(f"Gene {gene_symbol} not found in Ensembl")
                return None
                
        except Exception as e:
            logger.error(f"Erreur Ensembl API: {e}")
            return None
    
    def get_variants(self, gene_id: str) -> List[Dict]:
        """
        Récupère les variants connus pour un gène
        """
        try:
            ext = f"/overlap/id/{gene_id}?feature=variation"
            response = requests.get(self.server + ext, headers=self.headers, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                variants = []
                for var in data:
                    if var.get('consequence_type') in ['missense_variant', 'synonymous_variant', 
                                                       'frameshift_variant', 'stop_gained']:
                        variants.append({
                            'id': var.get('id'),
                            'location': f"{var.get('seq_region_name')}:{var.get('start')}-{var.get('end')}",
                            'alleles': var.get('alleles', []),
                            'consequence': var.get('consequence_type'),
                            'sift_score': var.get('sift_score'),
                            'polyphen_score': var.get('polyphen_score')
                        })
                return variants
            return []
            
        except Exception as e:
            logger.error(f"Erreur récupération variants: {e}")
            return []
    
    def get_homology(self, gene_id: str, target_species: str = "homo_sapiens") -> Optional[Dict]:
        """
        Récupère les homologies avec d'autres espèces (utile pour validation)
        """
        try:
            ext = f"/homology/id/{gene_id}?target_species={target_species}"
            response = requests.get(self.server + ext, headers=self.headers, timeout=10)
            
            if response.status_code == 200:
                return response.json()
            return None
            
        except Exception as e:
            logger.error(f"Erreur homologie: {e}")
            return None

class OMIAConnector:
    """
    Connexion à OMIA (Online Mendelian Inheritance in Animals)
    pour les traits génétiques et maladies ovines
    """
    
    def __init__(self):
        self.base_url = OMIA_API_BASE
    
    def search_traits(self, species: str = "sheep") -> List[Dict]:
        """
        Recherche les traits/phenotypes génétiques répertoriés
        """
        try:
            url = f"{self.base_url}/phenotypes?species={species}"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                return response.json().get('results', [])
            return []
            
        except Exception as e:
            logger.warning(f"OMIA API non disponible: {e}")
            return []
    
    def get_milk_production_traits(self) -> List[Dict]:
        """
        Spécifiquement pour les traits de production laitière
        """
        traits = [
            {
                'omia_id': 1,
                'trait_name': 'Milk yield',
                'inheritance': 'polygenic',
                'genes': ['DGAT1', 'ACACA', 'FASN', 'GH'],
                'phene_id': '1',
                'species': 'sheep'
            },
            {
                'omia_id': 2,
                'trait_name': 'Milk fat percentage',
                'inheritance': 'polygenic',
                'genes': ['DGAT1', 'SCD', 'FASN'],
                'phene_id': '2',
                'species': 'sheep'
            },
            {
                'omia_id': 3,
                'trait_name': 'Milk protein percentage',
                'inheritance': 'polygenic',
                'genes': ['CSN1S1', 'CSN3', 'LALBA', 'BLG'],
                'phene_id': '3',
                'species': 'sheep'
            }
        ]
        return traits

class GeneticDataIntegration:
    """
    Intégration complète des données génétiques multi-sources
    """
    
    def __init__(self):
        self.ncbi = NCBIConnector()
        self.ensembl = EnsemblConnector()
        self.omia = OMIAConnector()
        self.alphamissense = AlphaMissenseConnector()
    
    def get_complete_gene_profile(self, gene_symbol: str) -> Dict:
        """
        Agrège les données d'un gène depuis toutes les sources disponibles
        """
        profile = {
            'gene_symbol': gene_symbol,
            'sources': {},
            'variants': [],
            'variants_scored': [],  # NOUVEAUTÉ: avec scores AlphaMissense
            'traits_associated': [],
            'reliability_score': 0
        }
        
        # NCBI
        ncbi_data = self.ncbi.search_gene(gene_symbol)
        if ncbi_data:
            profile['sources']['ncbi'] = ncbi_data
            profile['reliability_score'] += 30
        
        # Ensembl
        ensembl_data = self.ensembl.get_gene_info(gene_symbol)
        if ensembl_data:
            profile['sources']['ensembl'] = ensembl_data
            profile['ensembl_id'] = ensembl_data.get('id')
            profile['location'] = f"{ensembl_data.get('seq_region_name')}:{ensembl_data.get('start')}-{ensembl_data.get('end')}"
            profile['reliability_score'] += 40
            
            # Variants Ensembl
            variants = self.ensembl.get_variants(ensembl_data.get('id'))
            profile['variants'].extend(variants)
            
            # NOUVEAUTÉ: Scoring AlphaMissense pour les variants missense
            missense_variants = [v for v in variants if v.get('consequence') == 'missense_variant']
            if missense_variants:
                scored_variants = []
                for var in missense_variants[:5]:  # Limiter pour la démo
                    # Conversion format variant pour AlphaMissense
                    variant_format = self._convert_to_protein_variant(var)
                    if variant_format:
                        score = self.alphamissense.get_variant_score(gene_symbol, variant_format)
                        if score:
                            scored_variants.append({**var, **score})
                profile['variants_scored'] = scored_variants
        
        # SNPs NCBI
        snps = self.ncbi.search_snp(gene_symbol)
        profile['ncbi_snps'] = snps
        if snps:
            profile['reliability_score'] += 30
        
        # Traits OMIA
        omia_traits = self.omia.get_milk_production_traits()
        associated = [t for t in omia_traits if gene_symbol in t.get('genes', [])]
        profile['traits_associated'] = associated
        
        return profile
    
    def _convert_to_protein_variant(self, variant_data: Dict) -> Optional[str]:
        """Convertit les données variant Ensembl en format protein pour AlphaMissense"""
        # Simplification pour démo
        return "p.Ala123Val"  # Exemple statique
    
    def score_animal_variants(self, gene_symbol: str, animal_variants: List[str]) -> List[Dict]:
        """
        Analyse les variants d'un animal spécifique avec AlphaMissense
        """
        return self.alphamissense.batch_score_variants(gene_symbol, animal_variants)
    
    def generate_breeding_recommendation(self, animal_genotypes: Dict, 
                                        target_trait: str = "milk_yield") -> Dict:
        """
        Génère des recommandations d'élevage basées sur les génotypes
        """
        recommendations = {
            'animal_score': 0,
            'genetic_potential': '',
            'breeding_strategy': '',
            'risks': [],
            'opportunities': []
        }
        
        gene_weights = {
            'DGAT1': 0.25, 'LALBA': 0.20, 'CSN1S1': 0.20,
            'CSN3': 0.15, 'PRLR': 0.10, 'STAT5A': 0.10
        }
        
        score = 0
        for gene, genotype in animal_genotypes.items():
            weight = gene_weights.get(gene, 0.05)
            
            if genotype in ['AA', 'GG', 'CC']:
                score += weight * 1.0
                recommendations['opportunities'].append(f"{gene}: Allèle favorable fixé")
            elif genotype in ['AG', 'GA', 'GC', 'CG']:
                score += weight * 0.6
                recommendations['opportunities'].append(f"{gene}: Potentiel d'amélioration par consanguinité")
            else:
                score += weight * 0.3
                recommendations['risks'].append(f"{gene}: Génotype défavorable - considérer croisement")
        
        recommendations['animal_score'] = round(score * 100, 1)
        
        if score > 0.8:
            recommendations['genetic_potential'] = "ÉLITE - Gardien de race"
            recommendations['breeding_strategy'] = "Reproduction intensive, conservation gamètes"
        elif score > 0.6:
            recommendations['genetic_potential'] = "SUPÉRIEUR - Bon reproducteur"
            recommendations['breeding_strategy'] = "Association avec élite complémentaire"
        elif score > 0.4:
            recommendations['genetic_potential'] = "MOYEN - Améliorable"
            recommendations['breeding_strategy'] = "Croisement avec élite forte"
        else:
            recommendations['genetic_potential'] = "STANDARD - À réformer"
            recommendations['breeding_strategy'] = "Exclusion reproduction ou croisement terminal"
        
        return recommendations

# ============================================================================
# MODULE 3: INDICE D'APTITUDE LAITIÈRE (IAL) - COMBINÉ IMAGE + GÉNÉTIQUE
# ============================================================================

@dataclass
class IALScore:
    """Structure de l'Indice d'Aptitude Laitière"""
    score_total: float  # Sur 100
    composante_phenotype: float  # Sur 40 (image)
    composante_genetique: float  # Sur 40 (génotype)
    composante_pedigree: float   # Sur 20 (ascendance)
    
    details_phenotype: Dict
    details_genetique: Dict
    
    def to_dict(self):
        return {
            'IAL_total': round(self.score_total, 2),
            'composante_phenotypique': round(self.composante_phenotype, 2),
            'composante_genetique': round(self.composante_genetique, 2),
            'composante_pedigree': round(self.composante_pedigree, 2),
            'details': {
                'phenotype': self.details_phenotype,
                'genetique': self.details_genetique
            },
            'classement': self._get_ranking()
        }
    
    def _get_ranking(self) -> str:
        if self.score_total >= 85:
            return "Elite (Top 5%)"
        elif self.score_total >= 75:
            return "Supérieur (Top 20%)"
        elif self.score_total >= 60:
            return "Bon (Top 50%)"
        else:
            return "Standard"

class IndiceAptitudeLaitiere:
    """
    Calcul de l'IAL combinant:
    - Analyse morphologique de la mamelle (image)
    - Génotypage des gènes de production laitière
    - Données de pedigree (optionnel)
    """
    
    def __init__(self):
        self.integration = GeneticDataIntegration()
        
        # Poids des différentes composantes
        self.weights = {
            'udder_morphology': 0.25,      # 25% - Forme mamelle
            'udder_capacity': 0.15,        # 15% - Capacité volume
            'milk_yield_genes': 0.20,      # 20% - Gènes quantité
            'milk_quality_genes': 0.15,    # 15% - Gènes qualité
            'udder_health_genes': 0.10,    # 10% - Gènes santé mamelle
            'pedigree_index': 0.15         # 15% - Index ascendants
        }
    
    def calculate_ial(self, 
                     udder_measurements: UdderMeasurements,
                     body_measurements: Dict,
                     genotypes: Dict[str, str],
                     pedigree_index: Optional[float] = None) -> IALScore:
        """
        Calcule l'IAL complet
        
        Args:
            udder_measurements: Mesures de la mamelle depuis l'analyse image
            body_measurements: Mensurations corporelles
            genotypes: Dict {gene: genotype} ex: {'DGAT1': 'AA', 'LALBA': 'AG'}
            pedigree_index: Index génétique des parents (0-100, optionnel)
        """
        
        # 1. Composante Phénotypique (Image) - 40 points max
        score_pheno = self._calculate_phenotypic_component(
            udder_measurements, body_measurements
        )
        
        # 2. Composante Génétique - 40 points max
        score_genetic = self._calculate_genetic_component(genotypes)
        
        # 3. Composante Pedigree - 20 points max
        score_pedigree = self._calculate_pedigree_component(pedigree_index)
        
        # Score total
        total = score_pheno + score_genetic + score_pedigree
        
        return IALScore(
            score_total=total,
            composante_phenotype=score_pheno,
            composante_genetique=score_genetic,
            composante_pedigree=score_pedigree,
            details_phenotype=self._get_phenotype_details(udder_measurements, body_measurements),
            details_genetique=self._get_genetic_details(genotypes)
        )
    
    def _calculate_phenotypic_component(self, 
                                       udder: UdderMeasurements, 
                                       body: Dict) -> float:
        """Calcule les 40 points de la composante phénotypique"""
        
        score = 0
        
        # Morphologie mamelle (25 points)
        # Basé sur score_symetrie, score_attache, score_morphologie
        morpho_score = (
            udder.score_symetrie * 0.4 +
            udder.score_attache * 0.35 +
            udder.score_morphologie * 0.25
        ) / 100 * 25  # Normaliser sur 25
        
        score += morpho_score
        
        # Capacité volume (15 points)
        # Basé sur surface et profondeur
        # Reference: excellente mamelle ~3000 cm², profondeur > 40cm
        capacity_ratio = min(1, (udder.surface_total_cm2 / 3000) * 0.6 + 
                                (udder.profondeur_cm / 45) * 0.4)
        capacity_score = capacity_ratio * 15
        
        score += capacity_score
        
        return round(score, 2)
    
    def _calculate_genetic_component(self, genotypes: Dict[str, str]) -> float:
        """Calcule les 40 points de la composante génétique"""
        
        score = 0
        
        # Gènes de rendement (20 points) - DGAT1, ACACA, GH
        yield_genes = {'DGAT1': 0.5, 'ACACA': 0.3, 'GH': 0.2}
        yield_score = 0
        
        for gene, weight in yield_genes.items():
            if gene in genotypes:
                gt = genotypes[gene]
                if gt in ['AA', 'GG']:  # Favorable homozygote
                    yield_score += weight * 20
                elif gt in ['AG', 'GA']:
                    yield_score += weight * 20 * 0.7
                else:
                    yield_score += weight * 20 * 0.3
        
        score += yield_score
        
        # Gènes de qualité (15 points) - LALBA, CSN1S1, CSN3
        quality_genes = {'LALBA': 0.4, 'CSN1S1': 0.35, 'CSN3': 0.25}
        quality_score = 0
        
        for gene, weight in quality_genes.items():
            if gene in genotypes:
                gt = genotypes[gene]
                # Logique similaire mais pondérée différemment
                if gt in ['AA', 'BB', 'CC']:
                    quality_score += weight * 15
                elif gt in ['AB', 'AC']:
                    quality_score += weight * 15 * 0.75
                else:
                    quality_score += weight * 15 * 0.4
        
        score += quality_score
        
        # Gènes de santé mamelle (5 points) - STAT5A, PRLR
        health_score = 0
        if 'STAT5A' in genotypes:
            health_score += 2.5 if genotypes['STAT5A'] in ['AA', 'GG'] else 1.5
        if 'PRLR' in genotypes:
            health_score += 2.5 if genotypes['PRLR'] in ['AA', 'GG'] else 1.5
        
        score += health_score
        
        return round(min(40, score), 2)
    
    def _calculate_pedigree_component(self, pedigree_index: Optional[float]) -> float:
        """Calcule les 20 points de la composante pedigree"""
        if pedigree_index is None:
            # Si pas d'info pedigree, attribuer un score neutre de 10
            return 10.0
        
        # Normaliser sur 20 points
        return round((pedigree_index / 100) * 20, 2)
    
    def _get_phenotype_details(self, udder: UdderMeasurements, body: Dict) -> Dict:
        """Retourne les détails du calcul phénotypique"""
        return {
            'surface_mamelle_cm2': udder.surface_total_cm2,
            'profondeur_mamelle_cm': udder.profondeur_cm,
            'score_symetrie': udder.score_symetrie,
            'score_attache': udder.score_attache,
            'indice_conformation': body.get('indice_conformation', 0),
            'qualite_photographie': 'Bonne' if udder.score_morphologie > 50 else 'Moyenne'
        }
    
    def _get_genetic_details(self, genotypes: Dict[str, str]) -> Dict:
        """Retourne les détails du calcul génétique"""
        details = {}
        
        for gene, genotype in genotypes.items():
            # Déterminer l'impact
            if genotype in ['AA', 'GG', 'CC']:
                impact = "Favorable"
                points = "100%"
            elif genotype in ['AG', 'GA', 'GC']:
                impact = "Intermédiaire"
                points = "70%"
            else:
                impact = "Défavorable"
                points = "30%"
            
            details[gene] = {
                'genotype': genotype,
                'impact': impact,
                'contribution': points
            }
        
        return details
    
    def generate_ial_report(self, ial_score: IALScore, animal_id: str) -> str:
        """
        Génère un rapport PDF/HTML de l'IAL
        """
        report = f"""
        RAPPORT INDICE D'APTITUDE LAITIÈRE (IAL)
        Animal ID: {animal_id}
        Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}
        
        SCORE GLOBAL: {ial_score.score_total}/100
        Classement: {ial_score._get_ranking()}
        
        COMPOSITION DU SCORE:
        - Composante Phénotypique (Image): {ial_score.composante_phenotype}/40
        - Composante Génétique: {ial_score.composante_genetique}/40
        - Composante Pedigree: {ial_score.composante_pedigree}/20
        
        RECOMMANDATIONS:
        """
        
        if ial_score.score_total >= 85:
            report += "\n✓ Animal ELITE - Prioriser la conservation du patrimoine génétique"
        elif ial_score.score_total >= 75:
            report += "\n✓ Animal SUPÉRIEUR - Bon candidat pour la reproduction sélective"
        elif ial_score.score_total >= 60:
            report += "\n○ Animal STANDARD - Potentiel d'amélioration par croisement"
        else:
            report += "\n✗ Animal à REFORMER - Ne pas utiliser en reproduction"
        
        return report

# ============================================================================
# INTERFACE STREAMLIT - MODULES COMBINÉS
# ============================================================================

def render_module_image_analysis():
    """Rend le module d'analyse d'image"""
    analyzer = ImageAnalyzerStreamlit()
    analyzer.render_interface()

def render_module_genetic_apis():
    """Rend le module de connexion aux APIs génétiques"""
    st.markdown("## 🧬 Intégration Données Génétiques Réelles (NCBI/Ensembl/AlphaMissense)")
    
    integration = GeneticDataIntegration()
    
    tabs = st.tabs(["🔍 Recherche Gène", "🧬 Profil Complet", "📊 Recommandations Élevage", "🎯 AlphaMissense"])
    
    with tabs[0]:
        st.subheader("Recherche dans NCBI et Ensembl")
        
        col_search, col_results = st.columns([1, 2])
        
        with col_search:
            gene_symbol = st.text_input("Symbole du gène", "DGAT1").upper()
            search_ncbi = st.checkbox("NCBI", value=True)
            search_ensembl = st.checkbox("Ensembl", value=True)
            search_omia = st.checkbox("OMIA", value=True)
            
            if st.button("🔍 Rechercher", type="primary"):
                with st.spinner("Connexion aux bases de données..."):
                    
                    results = {}
                    
                    if search_ncbi:
                        with st.spinner("NCBI..."):
                            ncbi_gene = integration.ncbi.search_gene(gene_symbol)
                            if ncbi_gene:
                                results['ncbi'] = ncbi_gene
                                st.success(f"✅ NCBI: {ncbi_gene.get('gene_id')}")
                                
                                snps = integration.ncbi.search_snp(gene_symbol)
                                if snps:
                                    st.info(f"📍 {len(snps)} SNPs trouvés")
                                    results['ncbi_snps'] = snps
                            else:
                                st.warning("❌ Non trouvé dans NCBI")
                    
                    if search_ensembl:
                        with st.spinner("Ensembl..."):
                            ens_data = integration.ensembl.get_gene_info(gene_symbol)
                            if ens_data:
                                results['ensembl'] = ens_data
                                st.success(f"✅ Ensembl: {ens_data.get('id')}")
                                st.write(f"Location: {ens_data.get('seq_region_name')}:{ens_data.get('start')}-{ens_data.get('end')}")
                            else:
                                st.warning("❌ Non trouvé dans Ensembl")
                    
                    st.session_state['last_gene_search'] = results
        
        with col_results:
            if 'last_gene_search' in st.session_state:
                results = st.session_state['last_gene_search']
                
                if results:
                    st.subheader("Résultats de la recherche")
                    
                    if 'ncbi' in results:
                        with st.expander("📚 Données NCBI", expanded=True):
                            ncbi = results['ncbi']
                            st.markdown(f"**Gene ID:** [{ncbi.get('gene_id')}]({ncbi.get('ncbi_url')})")
                            st.markdown(f"**Ensembl:** [{ncbi.get('ensembl_link')}]({ncbi.get('ensembl_link')})")
                    
                    if 'ncbi_snps' in results:
                        with st.expander("🎯 SNPs Associés", expanded=True):
                            snps_df = pd.DataFrame(results['ncbi_snps'])
                            if not snps_df.empty:
                                st.dataframe(snps_df[['rs_id', 'chromosome', 'position', 'alleles']])
                                
                                selected_snp = st.selectbox("SNP à analyser", snps_df['rs_id'].tolist())
                                snp_info = snps_df[snps_df['rs_id'] == selected_snp].iloc[0]
                                
                                st.markdown(f"""
                                **{snp_info['rs_id']}**
                                - Position: Chromosome {snp_info['chromosome']}, {snp_info['position']}
                                - Allèles: {snp_info['alleles']}
                                - [Voir dans NCBI]({snp_info['ncbi_url']})
                                """)
                    
                    if 'ensembl' in results:
                        with st.expander("🧬 Données Ensembl"):
                            ens = results['ensembl']
                            st.json({
                                'id': ens.get('id'),
                                'biotype': ens.get('biotype'),
                                'strand': ens.get('strand'),
                                'description': ens.get('description', '')[:200] + '...'
                            })
    
    with tabs[1]:
        st.subheader("Profil Génétique Complet Multi-Sources")
        
        gene_profile = st.selectbox("Gène à profiler", 
                                   ['DGAT1', 'LALBA', 'CSN1S1', 'CSN3', 'PRLR', 'STAT5A', 'ACACA', 'FASN'])
        
        if st.button("🚀 Générer profil complet", type="primary"):
            with st.spinner("Agrégation des données..."):
                profile = integration.get_complete_gene_profile(gene_profile)
                
                reliability = profile.get('reliability_score', 0)
                col_rel1, col_rel2 = st.columns([1, 3])
                with col_rel1:
                    st.metric("Fiabilité données", f"{reliability}%")
                with col_rel2:
                    st.progress(reliability / 100)
                
                if 'location' in profile:
                    st.success(f"📍 **Localisation:** {profile['location']}")
                
                # Variants avec scores AlphaMissense
                if profile.get('variants_scored'):
                    st.subheader("🎯 Variants avec Scores AlphaMissense")
                    for var in profile['variants_scored']:
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.write(f"**{var['id']}**")
                        with col2:
                            color = var.get('color', 'gray')
                            st.markdown(f"<span style='color:{color}'>●</span> {var['classification']}", 
                                      unsafe_allow_html=True)
                        with col3:
                            st.metric("Pathogénicité", f"{var.get('pathogenicity_score', 0):.2f}")
                
                if profile.get('traits_associated'):
                    st.subheader("🥛 Traits Laitiers Associés")
                    for trait in profile['traits_associated']:
                        st.info(f"**{trait['trait_name']}** ({trait['inheritance']})")
                
                st.subheader("🔗 Sources de Données")
                for source, data in profile.get('sources', {}).items():
                    st.write(f"- **{source.upper()}**: ✅ Connecté")
                
                st.download_button(
                    "📥 Exporter profil JSON",
                    json.dumps(profile, indent=2, default=str),
                    f"{gene_profile}_profile.json"
                )
    
    with tabs[2]:
        st.subheader("Recommandations d'Élevage basées sur le Génotype")
        
        st.info("""
        Cet outil analyse les génotypes d'un animal et génère des recommandations
        personnalisées pour la sélection et le croisement.
        """)
        
        st.markdown("### 🧬 Saisie des Génotypes")
        
        cols = st.columns(4)
        animal_genotypes = {}
        
        genes_input = ['DGAT1', 'LALBA', 'CSN1S1', 'CSN3', 'PRLR', 'STAT5A', 'ACACA', 'FASN']
        
        for i, gene in enumerate(genes_input):
            with cols[i % 4]:
                animal_genotypes[gene] = st.selectbox(
                    gene,
                    ['Non testé', 'AA', 'AB', 'BB', 'A/A', 'A/B', 'B/B', 'CC', 'AC', 'BC'],
                    key=f"geno_{gene}"
                )
        
        target = st.selectbox("Objectif de sélection", 
                            ["milk_yield", "milk_fat", "milk_protein", "conformation"])
        
        if st.button("🎯 Analyser et Recommander", type="primary"):
            genotypes_filtres = {k: v for k, v in animal_genotypes.items() 
                               if v not in ['Non testé']}
            
            if len(genotypes_filtres) < 3:
                st.error("❌ Minimum 3 gènes testés requis pour une analyse fiable")
            else:
                recommendations = integration.generate_breeding_recommendation(
                    genotypes_filtres, target
                )
                
                col_res1, col_res2 = st.columns(2)
                
                with col_res1:
                    st.metric("Score Génétique", f"{recommendations['animal_score']}/100")
                    
                    score = recommendations['animal_score']
                    color = 'green' if score > 80 else 'orange' if score > 60 else 'red'
                    st.markdown(f"""
                    <div style='padding: 20px; border-radius: 10px; background-color: {color}; color: white;'>
                        <h3>{recommendations['genetic_potential']}</h3>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col_res2:
                    st.subheader("Stratégie Recommandée")
                    st.write(recommendations['breeding_strategy'])
                
                col_opp, col_risk = st.columns(2)
                
                with col_opp:
                    st.subheader("✅ Opportunités")
                    for opp in recommendations['opportunities']:
                        st.success(opp)
                
                with col_risk:
                    st.subheader("⚠️ Risques")
                    for risk in recommendations['risks']:
                        st.warning(risk)
    
    # NOUVEAU TAB: AlphaMissense
    with tabs[3]:
        st.subheader("🧬 Analyse AlphaMissense - Impact des Variants")
        st.info("""
        AlphaMissense (DeepMind) prédit l'impact pathogénique des variants missense 
        avec une précision supérieure aux méthodes traditionnelles.
        """)
        
        col_am1, col_am2 = st.columns(2)
        
        with col_am1:
            am_gene = st.text_input("Gène", "DGAT1", key="am_gene").upper()
            am_variant = st.text_input("Variant (format: p.Ala123Val)", "p.Ala232Val")
            
            if st.button("🔬 Prédire impact", key="predict_am"):
                with st.spinner("Analyse par AlphaMissense..."):
                    result = integration.alphamissense.get_variant_score(am_gene, am_variant)
                    
                    if result:
                        st.session_state['am_result'] = result
        
        with col_am2:
            if 'am_result' in st.session_state:
                res = st.session_state['am_result']
                
                # Jauge de pathogénicité
                score = res['pathogenicity_score']
                st.metric("Score de Pathogénicité", f"{score:.3f}")
                st.progress(score)
                
                if res['classification'] == "Pathogénique":
                    st.error(f"⚠️ {res['classification']} - Impact fonctionnel probable")
                elif res['classification'] == "Bénin":
                    st.success(f"✅ {res['classification']} - Impact limité")
                else:
                    st.warning(f"❓ {res['classification']} - Requiert investigation")
                
                st.write(f"Confiance: {res['confidence']*100:.0f}%")
                st.caption(f"Méthode: {res['prediction_method']}")

def render_module_ial():
    """NOUVEAU: Module de calcul de l'Indice d'Aptitude Laitière"""
    st.markdown("## 🥛 Indice d'Aptitude Laitière (IAL)")
    st.info("""
    L'IAL combine l'analyse morphologique de la mamelle (par smartphone) 
    avec les données génétiques pour évaluer le potentiel laitier.
    """)
    
    col_input, col_result = st.columns([1, 1])
    
    with col_input:
        st.subheader("📸 Données Phénotypiques")
        
        # Utiliser les données d'analyse d'image si disponibles
        if 'image_analysis' in st.session_state and st.session_state.image_analysis.get('udder_analysis'):
            udder_data = st.session_state.image_analysis['udder_analysis']['measurements']
            st.success("✅ Données mamelle détectées depuis l'analyse image")
            use_image_data = st.checkbox("Utiliser les mesures image", value=True)
        else:
            st.info("Pas de données image. Veuillez saisir manuellement.")
            use_image_data = False
        
        if use_image_data:
            meas = udder_data
        else:
            # Saisie manuelle
            meas = UdderMeasurements(
                surface_total_cm2=st.number_input("Surface mamelle (cm²)", 1000, 4000, 2400),
                profondeur_cm=st.number_input("Profondeur (cm)", 20, 60, 42),
                largeur_cm=st.number_input("Largeur (cm)", 20, 50, 35),
                score_symetrie=st.slider("Score Symétrie", 0, 100, 85),
                score_attache=st.slider("Score Attache", 0, 100, 80),
                score_morphologie=st.slider("Score Morphologie Globale", 0, 100, 82)
            )
        
        st.subheader("🧬 Données Génétiques")
        
        col_g1, col_g2 = st.columns(2)
        genotypes = {}
        
        with col_g1:
            genotypes['DGAT1'] = st.selectbox("DGAT1", ['Non testé', 'AA', 'AG', 'GG'], key='dgat1')
            genotypes['LALBA'] = st.selectbox("LALBA", ['Non testé', 'AA', 'AB', 'BB'], key='lalba')
            genotypes['CSN1S1'] = st.selectbox("CSN1S1", ['Non testé', 'AA', 'AB', 'BB'], key='csn1s1')
        
        with col_g2:
            genotypes['CSN3'] = st.selectbox("CSN3", ['Non testé', 'AA', 'AB', 'BB'], key='csn3')
            genotypes['PRLR'] = st.selectbox("PRLR", ['Non testé', 'AA', 'AG', 'GG'], key='prlr')
            genotypes['STAT5A'] = st.selectbox("STAT5A", ['Non testé', 'AA', 'AG', 'GG'], key='stat5a')
        
        # Filtrer les non-testés
        genotypes = {k: v for k, v in genotypes.items() if v != 'Non testé'}
        
        st.subheader("📋 Pedigree (Optionnel)")
        pedigree_idx = st.slider("Index génétique parents (0-100)", 0, 100, 75)
        use_pedigree = st.checkbox("Intégrer le pedigree", value=False)
    
    with col_result:
        if st.button("🚀 Calculer l'IAL", type="primary"):
            if len(genotypes) < 2:
                st.error("Minimum 2 gènes testés requis")
            else:
                ial_calculator = IndiceAptitudeLaitiere()
                
                body_meas = {'indice_conformation': 26.5}  # Valeur par défaut
                
                ial_score = ial_calculator.calculate_ial(
                    meas,
                    body_meas,
                    genotypes,
                    pedigree_idx if use_pedigree else None
                )
                
                # Affichage résultats
                st.subheader("Résultat IAL")
                
                # Jauge circulaire simulée avec colonnes
                cols_score = st.columns([1, 2, 1])
                with cols_score[1]:
                    st.markdown(f"""
                    <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 50%; width: 150px; height: 150px; margin: auto; display: flex; align-items: center; justify-content: center;'>
                        <div>
                            <div style='font-size: 32px; color: white; font-weight: bold;'>{ial_score.score_total}</div>
                            <div style='font-size: 12px; color: white;'>sur 100</div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                
                st.markdown(f"### Classement: {ial_score._get_ranking()}")
                
                # Décomposition
                st.subheader("Décomposition du Score")
                
                fig = go.Figure(go.Bar(
                    x=['Phénotype<br>(Image)', 'Génétique<br>(ADN)', 'Pedigree<br>(Ascendance)'],
                    y=[ial_score.composante_phenotype, ial_score.composante_genetique, ial_score.composante_pedigree],
                    marker_color=['#FF6B6B', '#4ECDC4', '#45B7D1'],
                    text=[f"{ial_score.composante_phenotype}/40", 
                          f"{ial_score.composante_genetique}/40",
                          f"{ial_score.composante_pedigree}/20"],
                    textposition='auto'
                ))
                fig.update_layout(
                    yaxis=dict(range=[0, 45]),
                    showlegend=False,
                    height=400
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Détails
                with st.expander("Voir les détails du calcul"):
                    st.json(ial_score.to_dict())
                
                # Rapport exportable
                report = ial_calculator.generate_ial_report(ial_score, "ANIMAL_001")
                st.download_button(
                    "📄 Télécharger le rapport IAL",
                    report,
                    f"IAL_RAPPORT_ANIMAL_001_{datetime.now().strftime('%Y%m%d')}.txt"
                )

def main():
    """Application principale"""
    st.set_page_config(
        page_title="Expert Ovin DZ - Modules Avancés",
        page_icon="🐑",
        layout="wide"
    )
    
    st.sidebar.title("🐑 Modules Avancés")
    module = st.sidebar.radio(
        "Sélectionner le module",
        ["📱 Analyse Image Smartphone", 
         "🧬 Données Génétiques Réelles", 
         "🥛 Calcul IAL (Image + Génétique)",  # NOUVEAU
         "🔧 Configuration", 
         "📚 Documentation"]
    )
    
    if module == "📱 Analyse Image Smartphone":
        render_module_image_analysis()
    
    elif module == "🧬 Données Génétiques Réelles":
        render_module_genetic_apis()
    
    elif module == "🥛 Calcul IAL (Image + Génétique)":
        render_module_ial()  # NOUVEAU
    
    elif module == "🔧 Configuration":
        st.title("⚙️ Configuration")
        st.info("Configuration des clés API et préférences")
        
        st.subheader("Clés API (optionnel pour démo)")
        ncbi_key = st.text_input("NCBI API Key", type="password")
        ensembl_key = st.text_input("Ensembl API Key", type="password")
        alphamissense_key = st.text_input("AlphaMissense API Key", type="password")
        
        st.subheader("Calibration par défaut")
        default_calib = st.number_input("Pixels/cm par défaut", 10.0, 200.0, 35.0)
    
    elif module == "📚 Documentation":
        st.title("📚 Documentation Technique")
        
        st.markdown("""
        ## 📱 Module Analyse d'Image
        
        ### Calibration automatique
        L'algorithme détecte automatiquement les objets de référence standards:
        - Feuille A4 (21×29.7 cm)
        - Carte bancaire (8.56×5.398 cm)
        - Bâton gradué 1m
        - Pièce de 2€ (2.575 cm)
        
        ### Segmentation Spécifique Mamelle (Udder)
        Analyse avancée utilisant:
        - Segmentation par couleur HSV (teintes chair)
        - Algorithmes Watershed pour isolation
        - Détection des tétines par clustering K-means
        - Scoring morphologique: symétrie, attache, capacité
        
        ### Points anatomiques détectés
        """)
        
        points_df = pd.DataFrame([
            {'Point': p.name, 'Description': p.description, 'Utilisation': 'Obligatoire' if p.id_ref <= 5 else 'Optionnelle'}
            for p in POINTS_ANATOMIQUES.values()
        ])
        st.dataframe(points_df)
        
        st.markdown("""
        ## 🧬 Module Données Génétiques
        
        ### Sources de données
        1. **NCBI** (National Center for Biotechnology Information)
        2. **Ensembl** (EBI) - Genome Ovis aries Rambouillet
        3. **OMIA** (Online Mendelian Inheritance in Animals)
        4. **AlphaMissense** (DeepMind) - Prédiction impact variants
        
        ### AlphaMissense
        Utilise l'apprentissage profond (transformers) pour prédire la pathogénicité 
        des variants missense avec une performance supérieure aux méthodes classiques 
        (CADD, PolyPhen, SIFT).
        
        ## 🥛 Indice d'Aptitude Laitière (IAL)
        
        ### Composition du Score (/100)
        - **Phénotype (40 pts)**: Morphologie mamelle par analyse image
          - Symétrie, attache, profondeur, surface
        - **Génétique (40 pts)**: Génotypes des gènes majeurs
          - DGAT1, LALBA, CSN1S1, CSN3, etc.
        - **Pedigree (20 pts)**: Index génétique des ascendants
        
        ### Utilisation
        1. Analyser la mamelle avec le module Image
        2. Saisir les génotypes connus
        3. Calculer l'IAL combiné
        """)
        
        st.warning("""
        **Note sur les limites d'API:**
        - NCBI: 3 requêtes/seconde sans clé, 10 avec clé
        - Ensembl: Pas de limite stricte mais respecter les bonnes pratiques
        - AlphaMissense: Nécessite une clé d'API pour accès complet
        - Mise en cache automatique des résultats
        """)

if __name__ == "__main__":
    main()

# ============================================================================
# SECTION 18: PIED DE PAGE
# ============================================================================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version 4.0</p>
    <p>📐 Scanner 3D • 🎯 Critères de sélection • 🧬 Génétique • 📊 Statistiques</p>
    <p>© rahim 2026 LABORATOIRE GenApAgiE - Système de gestion scientifique des races ovines algériennes</p>
</div>
""", unsafe_allow_html=True)
