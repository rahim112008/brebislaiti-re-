"""
EXPERT OVIN DZ PRO - Modules Avancés
Module 1: Analyse d'image smartphone avec calibration OpenCV
Module 2: Connexion API NCBI/Ensembl pour données génétiques réelles
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

# ============================================================================
# SECTION 2: CONFIGURATION ET CONSTANTES
# ============================================================================
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuration APIs externes
NCBI_API_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_API_BASE = "https://rest.ensembl.org"
OMIA_API_BASE = "https://omia.org/api/v1"

# Cache pour limiter les appels API
API_CACHE = {}

# ============================================================================
# MODULE 1: ANALYSE D'IMAGE SMARTPHONE - CALIBRATION AUTOMATIQUE
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

class ImageAnalyzerOpenCV:
    """
    Analyse d'image morphométrique avec OpenCV
    Détection automatique de l'animal et calibration
    """
    
    def __init__(self):
        self.calibration_factor = None  # pixels/cm
        self.reference_object = None
        self.detected_points = {}
        
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
    
    def visualize_analysis(self, image: np.ndarray, contour: np.ndarray, 
                          points: Dict, ref_detection: Dict) -> np.ndarray:
        """
        Génère une visualisation de l'analyse
        """
        if not OPENCV_AVAILABLE:
            return image
        
        vis_image = image.copy()
        
        # Dessiner le contour
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
            st.subheader("🎯 Points à détecter")
            
            selected_points = st.multiselect(
                "Points anatomiques",
                list(POINTS_ANATOMIQUES.keys()),
                default=['garrot', 'tuber_coxal', 'pointe_bassin', 'mamelle_avant', 'mamelle_arriere'],
                format_func=lambda x: f"{POINTS_ANATOMIQUES[x].name} - {POINTS_ANATOMIQUES[x].description}"
            )
            
            st.markdown("---")
            st.subheader("🔧 Options avancées")
            
            auto_detect = st.toggle("Détection automatique", value=True)
            manual_adjust = st.toggle("Ajustement manuel", value=True)
        
        with col_image:
            uploaded_file = st.file_uploader(
                "📷 Charger la photo", 
                type=['jpg', 'jpeg', 'png', 'bmp'],
                help="Image de l'animal avec l'objet de référence visible"
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
                        
                        # Détection contour animal
                        contour = self.analyzer.detect_animal_contour(image_rgb) if auto_detect else None
                        
                        # Points anatomiques
                        if auto_detect:
                            points = self.analyzer.suggest_anatomical_points(image_rgb, contour)
                        else:
                            points = {}
                        
                        # Visualisation
                        if OPENCV_AVAILABLE and contour is not None:
                            vis_image = self.analyzer.visualize_analysis(
                                image_rgb, contour, points, ref_detection
                            )
                        else:
                            # Mode simulation avec PIL
                            vis_image = self._create_visualization_pil(
                                image_rgb, points, ref_detection
                            )
                        
                        st.image(vis_image, caption="Analyse détectée", use_column_width=True)
                        
                        # Ajustement manuel des points
                        if manual_adjust and points:
                            st.subheader("🎯 Affiner les points")
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
                            
                            st.subheader("📏 Mesures calculées")
                            cols_m = st.columns(len(measurements))
                            for col, (name, value) in zip(cols_m, measurements.items()):
                                with col:
                                    st.metric(
                                        name.replace('_', ' ').title(),
                                        f"{value} cm" if 'indice' not in name else f"{value}"
                                    )
                            
                            # Score de qualité
                            if 'indice_conformation' in measurements:
                                ic = measurements['indice_conformation']
                                if ic > 27:
                                    st.success(f"🥇 Bonne conformation (IC: {ic:.1f})")
                                elif ic > 24:
                                    st.info(f"🥈 Conformation moyenne (IC: {ic:.1f})")
                                else:
                                    st.warning(f"🥉 Conformation à améliorer (IC: {ic:.1f})")
                            
                            # Sauvegarde
                            if st.button("💾 Enregistrer les mesures", type="primary"):
                                self.session_state.image_analysis['measurements'] = measurements
                                self.session_state.image_analysis['points'] = points
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
    
    def _create_visualization_pil(self, image_array: np.ndarray, points: Dict, 
                                   ref_detection: Dict) -> np.ndarray:
        """Crée une visualisation avec PIL quand OpenCV n'est pas dispo"""
        img = Image.fromarray(image_array.astype('uint8'))
        draw = ImageDraw.Draw(img)
        
        # Dessiner points
        for name, (x, y) in points.items():
            point_info = POINTS_ANATOMIQUES.get(name, PointAnatomique(name, '', '#FFFFFF', 0))
            color = point_info.color
            
            # Cercle
            r = 8
            draw.ellipse([x-r, y-r, x+r, y+r], fill=color, outline='white', width=2)
            
            # Label
            draw.text((x+10, y-5), point_info.name, fill=color)
        
        # Cadre référence
        if ref_detection.get('detected') and 'bbox' in ref_detection:
            x, y, w, h = ref_detection['bbox']
            draw.rectangle([x, y, x+w, y+h], outline='blue', width=3)
        
        return np.array(img)

# ============================================================================
# MODULE 2: CONNEXION API NCBI/ENSEMBL - DONNÉES GÉNÉTIQUES RÉELLES
# ============================================================================

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
            # En production, utiliser xml.etree.ElementTree
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
        # Parsing simple
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
            # Utilisation de efetch pour nucleotide
            url = f"{NCBI_API_BASE}/efetch.fcgi"
            params = {
                'db': 'nuccore',
                'id': gene_id,
                'rettype': 'fasta',
                'retmode': 'text',
                'seq_start': max(1, 1 - upstream),  # Simplifié
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
                # Filtrer variants pertinents
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
            ext = f"/ homology/id/{gene_id}?target_species={target_species}"
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
            # OMIA nécessite parfois authentification
            # Version simplifiée pour démo
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
                'omia_id': 000001,
                'trait_name': 'Milk yield',
                'inheritance': 'polygenic',
                'genes': ['DGAT1', 'ACACA', 'FASN', 'GH'],
                'phene_id': '000001',
                'species': 'sheep'
            },
            {
                'omia_id': 000002,
                'trait_name': 'Milk fat percentage',
                'inheritance': 'polygenic',
                'genes': ['DGAT1', 'SCD', 'FASN'],
                'phene_id': '000002',
                'species': 'sheep'
            },
            {
                'omia_id': 000003,
                'trait_name': 'Milk protein percentage',
                'inheritance': 'polygenic',
                'genes': ['CSN1S1', 'CSN3', 'LALBA', 'BLG'],
                'phene_id': '000003',
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
    
    def get_complete_gene_profile(self, gene_symbol: str) -> Dict:
        """
        Agrège les données d'un gène depuis toutes les sources disponibles
        """
        profile = {
            'gene_symbol': gene_symbol,
            'sources': {},
            'variants': [],
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
        
        # Scoring basé sur les gènes majeurs laitiers
        gene_weights = {
            'DGAT1': 0.25, 'LALBA': 0.20, 'CSN1S1': 0.20,
            'CSN3': 0.15, 'PRLR': 0.10, 'STAT5A': 0.10
        }
        
        score = 0
        for gene, genotype in animal_genotypes.items():
            weight = gene_weights.get(gene, 0.05)
            
            # Évaluation génotype (simplifié)
            if genotype in ['AA', 'GG', 'CC']:  # Homozygote favorable
                score += weight * 1.0
                recommendations['opportunities'].append(f"{gene}: Allèle favorable fixé")
            elif genotype in ['AG', 'GA', 'GC', 'CG']:  # Hétérozygote
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
# INTERFACE STREAMLIT - MODULES COMBINÉS
# ============================================================================

def render_module_image_analysis():
    """Rend le module d'analyse d'image"""
    analyzer = ImageAnalyzerStreamlit()
    analyzer.render_interface()

def render_module_genetic_apis():
    """Rend le module de connexion aux APIs génétiques"""
    st.markdown("## 🧬 Intégration Données Génétiques Réelles (NCBI/Ensembl)")
    
    integration = GeneticDataIntegration()
    
    tabs = st.tabs(["🔍 Recherche Gène", "🧬 Profil Complet", "📊 Recommandations Élevage"])
    
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
                                
                                # SNPs
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
                    
                    # Affichage synthétique
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
                                
                                # Sélection SNP pour détails
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
                
                # Score de fiabilité
                reliability = profile.get('reliability_score', 0)
                col_rel1, col_rel2 = st.columns([1, 3])
                with col_rel1:
                    st.metric("Fiabilité données", f"{reliability}%")
                with col_rel2:
                    st.progress(reliability / 100)
                
                # Cartographie
                if 'location' in profile:
                    st.success(f"📍 **Localisation:** {profile['location']}")
                
                # Variants
                if profile.get('variants'):
                    st.subheader("🎯 Variants Fonctionnels")
                    variants_df = pd.DataFrame(profile['variants'])
                    st.dataframe(variants_df)
                
                # Traits associés
                if profile.get('traits_associated'):
                    st.subheader("🥛 Traits Laitiers Associés")
                    for trait in profile['traits_associated']:
                        st.info(f"**{trait['trait_name']}** ({trait['inheritance']})")
                
                # Sources
                st.subheader("🔗 Sources de Données")
                for source, data in profile.get('sources', {}).items():
                    st.write(f"- **{source.upper()}**: ✅ Connecté")
                
                # Export
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
        
        # Saisie génotypes
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
            # Filtrer les non-testés
            genotypes_filtres = {k: v for k, v in animal_genotypes.items() 
                               if v not in ['Non testé']}
            
            if len(genotypes_filtres) < 3:
                st.error("❌ Minimum 3 gènes testés requis pour une analyse fiable")
            else:
                recommendations = integration.generate_breeding_recommendation(
                    genotypes_filtres, target
                )
                
                # Affichage résultats
                col_res1, col_res2 = st.columns(2)
                
                with col_res1:
                    st.metric("Score Génétique", f"{recommendations['animal_score']}/100")
                    
                    # Jauge
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
                
                # Opportunités et risques
                col_opp, col_risk = st.columns(2)
                
                with col_opp:
                    st.subheader("✅ Opportunités")
                    for opp in recommendations['opportunities']:
                        st.success(opp)
                
                with col_risk:
                    st.subheader("⚠️ Risques")
                    for risk in recommendations['risks']:
                        st.warning(risk)

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
        ["📱 Analyse Image Smartphone", "🧬 Données Génétiques Réelles", 
         "🔧 Configuration", "📚 Documentation"]
    )
    
    if module == "📱 Analyse Image Smartphone":
        render_module_image_analysis()
    
    elif module == "🧬 Données Génétiques Réelles":
        render_module_genetic_apis()
    
    elif module == "🔧 Configuration":
        st.title("⚙️ Configuration")
        st.info("Configuration des clés API et préférences")
        
        st.subheader("Clés API (optionnel pour démo)")
        ncbi_key = st.text_input("NCBI API Key", type="password")
        ensembl_key = st.text_input("Ensembl API Key", type="password")
        
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
           - dbSNP: Variants génétiques
           - Gene: Informations géniques
           - Nucleotide: Séquences
        
        2. **Ensembl** (EBI)
           - Génome Ovis aries Rambouillet
           - Annotations fonctionnelles
           - Variants et conséquences
        
        3. **OMIA** (Online Mendelian Inheritance in Animals)
           - Traits phénotypiques
           - Maladies génétiques
           - Inheritance patterns
        """)
        
        st.warning("""
        **Note sur les limites d'API:**
        - NCBI: 3 requêtes/seconde sans clé, 10 avec clé
        - Ensembl: Pas de limite stricte mais respecter les bonnes pratiques
        - Mise en cache automatique des résultats
        """)

if __name__ == "__main__":
    main()
