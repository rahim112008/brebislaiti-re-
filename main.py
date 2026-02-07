"""
EXPERT OVIN DZ PRO - VERSION ENTERPRISE 2026.02.1
Module Scanner IA Morphométrique avec Calibration Étalon
Module principal de scan morphométrique par vision IA
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import logging
import hashlib
import io
import base64
from datetime import datetime, date, timedelta
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union, Literal, Any
from contextlib import contextmanager
from functools import lru_cache
from enum import Enum
import re

# Vision / Image processing
try:
    from PIL import Image, ImageDraw, ImageFont, ImageOps
    import cv2
    VISION_AVAILABLE = True
except ImportError:
    VISION_AVAILABLE = False
    logging.warning("OpenCV/Pillow non disponible - Mode simulation uniquement")

# Configuration logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler('ovin_pro.log'), logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION & CONSTANTES
# ============================================================================

class EtalonType(Enum):
    """Types d'étalons disponibles pour la calibration"""
    BATON_1M = "baton_1m"           # 100 cm
    CARTE_BANCAIRE = "carte_bancaire"  # 8.56 cm × 5.398 cm (ISO 7810 ID-1)
    FEUILLE_A4 = "feuille_a4"       # 29.7 cm × 21 cm

@dataclass
class AppConfig:
    """Configuration globale de l'application"""
    DB_PATH: str = "data/ovin_enterprise.db"
    CACHE_TTL: int = 3600
    ALIGNMENT_THRESHOLD: float = 85.0
    MIN_SEQ_LENGTH: int = 50
    
    # Dimensions étalons référence (cm)
    ETALON_DIMS: Dict[str, Dict[str, Any]] = field(default_factory=lambda: {
        EtalonType.BATON_1M.value: {
            "longueur_reelle_cm": 100.0,
            "largeur_cm": 2.5,  # bâton standard
            "nom": "Bâton de 1 mètre",
            "couleur_recommandee": "Jaune/Orange vif",
            "tolerance_cm": 0.5
        },
        EtalonType.CARTE_BANCAIRE.value: {
            "longueur_reelle_cm": 8.56,
            "largeur_reelle_cm": 5.398,
            "nom": "Carte bancaire (ISO 7810)",
            "couleur_recommandee": "Peu importe",
            "tolerance_cm": 0.1
        },
        EtalonType.FEUILLE_A4.value: {
            "longueur_reelle_cm": 29.7,
            "largeur_reelle_cm": 21.0,
            "nom": "Feuille A4 standard",
            "couleur_recommandee": "Blanche (contraste maximal)",
            "tolerance_cm": 0.3
        }
    })

CONFIG = AppConfig()

# ============================================================================
# DATABASE MANAGER (Manquant dans le code original)
# ============================================================================

class DatabaseManager:
    """Gestionnaire de base de données avec context manager"""
    
    def __init__(self, db_path: str = CONFIG.DB_PATH):
        self.db_path = db_path
        self.conn = None
        self._ensure_data_dir()
    
    def _ensure_data_dir(self):
        """Crée le répertoire data si inexistant"""
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
    
    def __enter__(self):
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def connect(self):
        """Établit la connexion à la base de données"""
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row
            logger.info(f"Connexion DB établie: {self.db_path}")
            return self.conn
        except sqlite3.Error as e:
            logger.error(f"Erreur connexion DB: {e}")
            raise
    
    def close(self):
        """Ferme la connexion à la base de données"""
        if self.conn:
            self.conn.close()
            self.conn = None
    
    def execute(self, query: str, params: Tuple = ()) -> sqlite3.Cursor:
        """Exécute une requête SQL"""
        if not self.conn:
            self.connect()
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            logger.error(f"Erreur SQL: {e}\nQuery: {query}")
            raise
    
    def fetch_one(self, query: str, params: Tuple = ()) -> Optional[sqlite3.Row]:
        """Récupère un seul résultat"""
        cursor = self.execute(query, params)
        return cursor.fetchone()
    
    def fetch_all(self, query: str, params: Tuple = ()) -> List[sqlite3.Row]:
        """Récupère tous les résultats"""
        cursor = self.execute(query, params)
        return cursor.fetchall()
    
    def fetch_df(self, query: str, params: Tuple = ()) -> pd.DataFrame:
        """Récupère les résultats dans un DataFrame pandas"""
        try:
            return pd.read_sql_query(query, self.conn, params=params)
        except Exception as e:
            logger.error(f"Erreur fetch_df: {e}")
            return pd.DataFrame()

# ============================================================================
# MOTEUR DE VISION & MORPHOMÉTRIE IA
# ============================================================================

class MorphoMetricsEngine:
    """
    Moteur de mesure morphométrique par analyse d'image
    Calibration par étalon + détection des points anatomiques clés
    """
    
    def __init__(self):
        self.reference_object_pixels = None
        self.pixels_per_cm = None
        self.etalon_type = None
        self.image_calibrated = None
        self.calibration_validated = False
        
    def set_etalon(self, etalon_type: EtalonType, pixel_length: float, 
                   pixel_width: Optional[float] = None) -> float:
        """
        Calibration: établit l'échelle pixels -> cm réels
        
        Args:
            etalon_type: Type d'étalon utilisé
            pixel_length: Longueur en pixels mesurée
            pixel_width: Largeur en pixels mesurée (optionnelle)
            
        Returns:
            Pixels par cm (ppm)
            
        Raises:
            ValueError: Si calibration impossible
        """
        self.etalon_type = etalon_type
        dims = CONFIG.ETALON_DIMS[etalon_type.value]
        
        try:
            # Calcul échelle basée sur la dimension la plus précise mesurée
            if pixel_width and etalon_type != EtalonType.BATON_1M:
                # Utilise les deux dimensions pour validation croisée
                ppm_long = pixel_length / dims["longueur_reelle_cm"]
                ppm_larg = pixel_width / dims["largeur_reelle_cm"]
                self.pixels_per_cm = (ppm_long + ppm_larg) / 2
                
                # Vérification cohérence (ratio doit matcher)
                ratio_mesure = pixel_length / pixel_width
                ratio_reel = dims["longueur_reelle_cm"] / dims["largeur_reelle_cm"]
                ecart_ratio = abs(ratio_mesure - ratio_reel) / ratio_reel
                
                if ecart_ratio > 0.1:  # 10% tolérance
                    logger.warning(f"Distorsion détectée: écart ratio {ecart_ratio:.1%}")
                    self.calibration_validated = False
                else:
                    self.calibration_validated = True
            else:
                # Pour bâton 1m, seule la longueur est utilisée
                self.pixels_per_cm = pixel_length / dims["longueur_reelle_cm"]
                self.calibration_validated = True
            
            self.reference_object_pixels = pixel_length
            logger.info(f"Calibration: {self.pixels_per_cm:.2f} px/cm avec {etalon_type.value}")
            
            return self.pixels_per_cm
            
        except ZeroDivisionError as e:
            logger.error(f"Erreur calibration: division par zéro - {e}")
            raise ValueError("Longueur d'étalon invalide (nulle)")
        except KeyError as e:
            logger.error(f"Erreur calibration: étalon inconnu - {e}")
            raise ValueError(f"Type d'étalon inconnu: {etalon_type}")
    
    def is_calibrated(self) -> bool:
        """Vérifie si le moteur est calibré"""
        return self.pixels_per_cm is not None and self.calibration_validated
    
    def pixels_to_cm(self, pixels: float) -> float:
        """Convertit une distance en pixels vers cm réels"""
        if not self.is_calibrated():
            raise ValueError("Calibration non effectuée ou invalide")
        return pixels / self.pixels_per_cm
    
    def cm_to_pixels(self, cm: float) -> float:
        """Convertit cm vers pixels (pour overlay)"""
        if not self.is_calibrated():
            raise ValueError("Calibration non effectuée ou invalide")
        return cm * self.pixels_per_cm
    
    def calculate_measurements(self, points: Dict[str, Tuple[int, int]]) -> Dict[str, float]:
        """
        Calcule les distances morphométriques à partir des points clics utilisateur
        
        Args:
            points: Dictionnaire des points anatomiques avec coordonnées (x, y)
            
        Returns:
            Dictionnaire des mesures calculées en cm
        """
        if not self.is_calibrated():
            return {"error": "Calibration requise"}
        
        measurements = {}
        
        try:
            # 1. LONGUEUR DU CORPS (Garrot -> Base de la queue)
            if all(k in points for k in ['garrot', 'base_queue']):
                dist_px = self._distance(points['garrot'], points['base_queue'])
                measurements['longueur_corps_cm'] = round(self.pixels_to_cm(dist_px), 1)
            
            # 2. HAUTEUR AU GARROT (Garrot -> Sol - approximé par patte)
            if all(k in points for k in ['garrot', 'patte_arriere']):
                # Hauteur = différence Y (image supposée de profil, niveau sol constant)
                hauteur_px = abs(points['garrot'][1] - points['patte_arriere'][1])
                measurements['hauteur_garrot_cm'] = round(self.pixels_to_cm(hauteur_px), 1)
            
            # 3. LARGEUR DU BASSIN (Hanche gauche -> Hanche droite - si vue arrière/disponible)
            hanche_keys = ['hanche_gauche', 'hanche_droite']
            if 'hanche' in points and len(hanche_keys) == 0:
                # Si seul 'hanche' est fourni (cas simple)
                pass
            elif all(k in points for k in hanche_keys):
                bassin_px = self._distance(points['hanche_gauche'], points['hanche_droite'])
                measurements['largeur_bassin_cm'] = round(self.pixels_to_cm(bassin_px), 1)
            
            # 4. CIRCONFÉRENCE DU CANON (Patte arrière - approximation elliptique)
            if 'canon' in points and 'canon_peripherie' in points:
                # Mesure diamètre, calcul circonférence = π × d
                diametre_px = self._distance(points['canon'], points['canon_peripherie']) * 2
                circonf_px = np.pi * diametre_px
                measurements['circonf_canon_cm'] = round(self.pixels_to_cm(circonf_px), 1)
            
            # 5. TOUR DE POITRINE (Poitrine -> Dos au niveau épaule)
            if all(k in points for k in ['poitrine', 'garrot']):
                # Approximation: mesure ligne droite, ajustement empirique ×1.2 pour circonférence
                poitrine_px = self._distance(points['poitrine'], points['garrot']) * 1.2
                measurements['tour_poitrine_cm'] = round(self.pixels_to_cm(poitrine_px), 1)
            
            # 6. MORPHOLOGIE MAMELLE
            mamelle_keys = ['mamelle_gauche', 'mamelle_droite', 'mamelle_arriere']
            if all(k in points for k in mamelle_keys):
                # Largeur attachment
                attache_px = self._distance(points['mamelle_gauche'], points['mamelle_droite'])
                measurements['attachment_mamelle_cm'] = round(self.pixels_to_cm(attache_px), 1)
                
                # Profondeur (arrière -> ligne attachment)
                center_x = (points['mamelle_gauche'][0] + points['mamelle_droite'][0]) // 2
                center_y = (points['mamelle_gauche'][1] + points['mamelle_droite'][1]) // 2
                profondeur_px = self._distance(points['mamelle_arriere'], (center_x, center_y))
                measurements['profondeur_mamelle_cm'] = round(self.pixels_to_cm(profondeur_px), 1)
                
                # Score calculé (formule simplifiée INRA)
                if 'attachment_mamelle_cm' in measurements and 'profondeur_mamelle_cm' in measurements:
                    score = min(10.0, max(1.0, 
                        (measurements['attachment_mamelle_cm'] / 10) + 
                        (measurements['profondeur_mamelle_cm'] / 5)
                    ))
                    measurements['score_mamelle'] = round(score, 1)
            
            # 7. INDICE COMPACTÉ (Poids estimé / Longueur)
            if 'longueur_corps_cm' in measurements and 'tour_poitrine_cm' in measurements:
                # Formule simplifiée: Poids estimé (kg) = (Tour poitrine² × Longueur) / 10800
                poids_estime = (measurements['tour_poitrine_cm'] ** 2 * measurements['longueur_corps_cm']) / 10800
                measurements['poids_estime_kg'] = round(poids_estime, 1)
                measurements['indice_compacite'] = round(poids_estime / measurements['longueur_corps_cm'], 2)
            
            return measurements
            
        except Exception as e:
            logger.error(f"Erreur calcul mesures: {e}")
            return {"error": f"Erreur calcul: {str(e)}"}
    
    def _distance(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> float:
        """Distance euclidienne entre deux points"""
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    def generate_overlay_image(self, image: Image.Image, points: Dict[str, Tuple[int, int]], 
                               measurements: Dict[str, float]) -> Image.Image:
        """Génère image annotée avec lignes de mesure"""
        if not VISION_AVAILABLE:
            raise ImportError("Pillow/OpenCV non disponible pour générer l'overlay")
        
        img_annotated = image.copy()
        draw = ImageDraw.Draw(img_annotated)
        
        # Couleurs standardisées
        COLOR_ETALON = (255, 0, 0)      # Rouge
        COLOR_MESURE = (0, 255, 0)      # Vert
        COLOR_POINT = (0, 0, 255)       # Bleu
        COLOR_TEXT = (255, 255, 255)    # Blanc
        COLOR_BG_TEXT = (0, 0, 0, 128)  # Noir semi-transparent
        
        try:
            font = ImageFont.truetype("arial.ttf", 20)
        except (OSError, AttributeError):
            font = ImageFont.load_default()
        
        # Dessine points cliqués
        for name, (x, y) in points.items():
            draw.ellipse([x-5, y-5, x+5, y+5], fill=COLOR_POINT, outline=COLOR_TEXT, width=2)
            # Fond pour texte lisible
            text_bbox = draw.textbbox((x+8, y-8), name, font=font)
            draw.rectangle(text_bbox, fill=COLOR_BG_TEXT)
            draw.text((x+8, y-8), name, fill=COLOR_TEXT, font=font)
        
        # Dessine lignes de mesure
        lines = [
            ('garrot', 'base_queue', 'longueur_corps_cm'),
            ('garrot', 'patte_arriere', 'hauteur_garrot_cm'),
            ('mamelle_gauche', 'mamelle_droite', 'attachment_mamelle_cm')
        ]
        
        for p1_name, p2_name, mesure_key in lines:
            if p1_name in points and p2_name in points:
                p1, p2 = points[p1_name], points[p2_name]
                draw.line([p1, p2], fill=COLOR_MESURE, width=3)
                
                # Label milieu avec fond
                mid_x = (p1[0] + p2[0]) // 2
                mid_y = (p1[1] + p2[1]) // 2
                if mesure_key in measurements:
                    text = f"{measurements[mesure_key]} cm"
                    text_bbox = draw.textbbox((mid_x, mid_y), text, font=font)
                    draw.rectangle(text_bbox, fill=COLOR_BG_TEXT)
                    draw.text((mid_x, mid_y), text, fill=COLOR_MESURE, font=font)
        
        return img_annotated

# ============================================================================
# INTERFACE SCANNER IA
# ============================================================================

def render_scanner_ia(db: DatabaseManager):
    """Interface complète de scan morphométrique avec calibration étalon"""
    st.header("📷 Scanner IA Morphométrique 1m")
    
    # Initialisation session state
    if 'etalon_type' not in st.session_state:
        st.session_state.etalon_type = None
    if 'morpho_engine' not in st.session_state:
        st.session_state.morpho_engine = MorphoMetricsEngine()
    if 'calibration_points' not in st.session_state:
        st.session_state.calibration_points = []
    if 'ppm' not in st.session_state:
        st.session_state.ppm = None
    
    st.markdown("""
    ### 📋 Protocole de Mesure Assistée par IA
    
    **Principe**: Calibration par objet de taille connue (étalon), puis mesures morphométriques 
    par détection des points anatomiques clés.
    
    **Précision attendue**: ±2 cm si protocole respecté
    """)
    
    # Étape 1: Choix de l'étalon
    st.subheader("1️⃣ Choix de l'Étalon de Calibration")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**🟨 Bâton 1 mètre** (Recommandé)")
        st.caption("Précision maximale. Placer à côté de l'animal, au sol.")
        if st.button("Sélectionner Bâton 1m", key="etalon_baton"):
            st.session_state.etalon_type = EtalonType.BATON_1M
    
    with col2:
        st.markdown("**💳 Carte Bancaire** (Urgence)")
        st.caption("ISO 7810: 8.56 × 5.40 cm. Placer sur le dos de l'animal.")
        if st.button("Sélectionner Carte", key="etalon_carte"):
            st.session_state.etalon_type = EtalonType.CARTE_BANCAIRE
    
    with col3:
        st.markdown("**📄 Feuille A4** (Standard)")
        st.caption("29.7 × 21 cm. Fixer verticalement près de l'animal.")
        if st.button("Sélectionner A4", key="etalon_a4"):
            st.session_state.etalon_type = EtalonType.FEUILLE_A4
    
    # Affichage étalon sélectionné
    if st.session_state.etalon_type:
        dims = CONFIG.ETALON_DIMS[st.session_state.etalon_type.value]
        st.success(f"✅ Étalonnage sélectionné: {dims['nom']} ({dims['longueur_reelle_cm']} cm)")
    
    # Étape 2: Upload image
    st.subheader("2️⃣ Acquisition de l'Image")
    
    uploaded_file = st.file_uploader(
        "📸 Photo de profil de l'animal avec l'étalon visible",
        type=['jpg', 'jpeg', 'png', 'bmp', 'tiff'],
        help="L'animal doit être de profil, tête à gauche, étalon visible et parallèle au corps"
    )
    
    if uploaded_file:
        try:
            # Lecture et traitement image
            image = Image.open(uploaded_file)
            
            # Conversion RGBA pour uniformité
            if image.mode != 'RGBA':
                image = image.convert('RGBA')
            
            # Redimensionnement pour affichage (conservation ratio)
            max_width = 800
            if image.width > max_width:
                ratio = max_width / image.width
                new_size = (max_width, int(image.height * ratio))
                display_image = image.resize(new_size, Image.Resampling.LANCZOS)
            else:
                display_image = image
            
            st.image(display_image, caption="Image source", use_container_width=True)
            
            if st.session_state.etalon_type:
                _render_calibration_step(display_image, db)
            else:
                st.warning("⚠️ Veuillez d'abord sélectionner un type d'étalon de calibration")
                
        except Exception as e:
            st.error(f"❌ Erreur lors du chargement de l'image: {str(e)}")
            logger.error(f"Erreur chargement image: {e}")
    
    # Section historique
    _render_history_section(db)

def _render_calibration_step(display_image: Image.Image, db: DatabaseManager):
    """Affiche l'interface de calibration"""
    st.subheader("3️⃣ Calibration sur l'Image")
    
    dims = CONFIG.ETALON_DIMS[st.session_state.etalon_type.value]
    st.markdown(f"""
    **Instructions**: 
    1. Cliquez sur **les deux extrémités** de l'étalon ({dims['nom']})
    2. Assurez-vous que la ligne est parallèle au sol
    """)
    
    # Interface de saisie calibration
    col_cal1, col_cal2 = st.columns([3, 1])
    
    with col_cal1:
        st.markdown("### 🎯 Coordonnées des extrémités de l'étalon")
        
        # Coordonnées étalon
        x1 = st.number_input("X1 (début étalon)", 0, display_image.width, 100, key="x1_cal")
        y1 = st.number_input("Y1 (début étalon)", 0, display_image.height, 400, key="y1_cal")
        x2 = st.number_input("X2 (fin étalon)", 0, display_image.width, 300, key="x2_cal")
        y2 = st.number_input("Y2 (fin étalon)", 0, display_image.height, 400, key="y2_cal")
        
        # Visualisation calibration
        calib_img = display_image.copy()
        draw_cal = ImageDraw.Draw(calib_img)
        draw_cal.line([(x1, y1), (x2, y2)], fill=(255, 0, 0), width=3)
        draw_cal.ellipse([x1-5, y1-5, x1+5, y1+5], fill=(255, 0, 0))
        draw_cal.ellipse([x2-5, y2-5, x2+5, y2+5], fill=(255, 0, 0))
        
        # Ajout texte
        font = ImageFont.load_default()
        draw_cal.text((x1+10, y1), "Début", fill=(255, 255, 255))
        draw_cal.text((x2+10, y2), "Fin", fill=(255, 255, 255))
        
        st.image(calib_img, caption="Vérification calibration", use_container_width=True)
    
    with col_cal2:
        pixel_distance = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        st.metric("Distance pixels", f"{pixel_distance:.1f} px")
        
        if st.button("✅ Valider Calibration", key="btn_validate_calib"):
            try:
                engine = st.session_state.morpho_engine
                ppm = engine.set_etalon(
                    st.session_state.etalon_type, 
                    pixel_distance
                )
                st.session_state.ppm = ppm
                st.success(f"Calibration: {ppm:.2f} pixels/cm")
                st.session_state.calibration_done = True
                st.rerun()
            except ValueError as e:
                st.error(f"Erreur calibration: {str(e)}")
            except Exception as e:
                st.error(f"Erreur inattendue: {str(e)}")
                logger.error(f"Erreur validation calibration: {e}")
    
    # Étape 4: Mesures morphométriques (si calibration validée)
    if getattr(st.session_state, 'calibration_done', False):
        _render_measurement_step(display_image, db)

def _render_measurement_step(display_image: Image.Image, db: DatabaseManager):
    """Affiche l'interface de mesure morphométrique"""
    st.subheader("4️⃣ Points de Mesure Morphométrique")
    
    st.markdown("""
    **Points à cliquer sur l'image**:
    1. **Garrot** (Point haut entre omoplates)
    2. **Base queue** (Insertion queue sur croupe)
    3. **Hanche** (Tuber coxal - pointe hanche)
    4. **Patte arrière** (Sol au niveau du jarret)
    5. **Canon** (Milieu du canon postérieur)
    6. **Poitrine** (Point le plus large de la poitrine)
    7. **Mamelles** (Attachment gauche/droite + arrière)
    """)
    
    # Interface de saisie des points
    points = {}
    
    with st.expander("📍 Saisie des Points Anatomiques", expanded=True):
        cols = st.columns(2)
        
        with cols[0]:
            st.markdown("**Points Gauche/Arrière**")
            points['garrot'] = (
                st.number_input("Garrot X", 0, display_image.width, 200, key="gx"),
                st.number_input("Garrot Y", 0, display_image.height, 150, key="gy")
            )
            points['base_queue'] = (
                st.number_input("Base Queue X", 0, display_image.width, 600, key="bqx"),
                st.number_input("Base Queue Y", 0, display_image.height, 200, key="bqy")
            )
            points['patte_arriere'] = (
                st.number_input("Patte Arrière X", 0, display_image.width, 550, key="pax"),
                st.number_input("Patte Arrière Y", 0, display_image.height, 500, key="pay")
            )
            points['mamelle_gauche'] = (
                st.number_input("Mamelle G X", 0, display_image.width, 300, key="mgx"),
                st.number_input("Mamelle G Y", 0, display_image.height, 450, key="mgy")
            )
        
        with cols[1]:
            st.markdown("**Points Droit/Devant**")
            points['hanche'] = (
                st.number_input("Hanche X", 0, display_image.width, 500, key="hx"),
                st.number_input("Hanche Y", 0, display_image.height, 250, key="hy")
            )
            points['poitrine'] = (
                st.number_input("Poitrine X", 0, display_image.width, 150, key="px"),
                st.number_input("Poitrine Y", 0, display_image.height, 350, key="py")
            )
            points['canon'] = (
                st.number_input("Canon X", 0, display_image.width, 550, key="cx"),
                st.number_input("Canon Y", 0, display_image.height, 480, key="cy")
            )
            points['mamelle_droite'] = (
                st.number_input("Mamelle D X", 0, display_image.width, 400, key="mdx"),
                st.number_input("Mamelle D Y", 0, display_image.height, 450, key="mdy")
            )
        
        # Points supplémentaires
        points['mamelle_arriere'] = (
            st.number_input("Mamelle Arrière X", 0, display_image.width, 350, key="max"),
            st.number_input("Mamelle Arrière Y", 0, display_image.height, 500, key="may")
        )
        points['canon_peripherie'] = (
            st.number_input("Périphérie Canon X", 0, display_image.width, 570, key="cpx"),
            st.number_input("Périphérie Canon Y", 0, display_image.height, 480, key="cpy")
        )
    
    if st.button("🧮 Calculer les Mesures", key="btn_calculate"):
        engine = st.session_state.morpho_engine
        measurements = engine.calculate_measurements(points)
        
        if "error" in measurements:
            st.error(measurements["error"])
        else:
            st.session_state.last_measurements = measurements
            st.session_state.last_points = points
            
            # Affichage résultats
            _display_results(measurements, points, display_image, db)

def _display_results(measurements: Dict[str, float], points: Dict[str, Tuple[int, int]], 
                     display_image: Image.Image, db: DatabaseManager):
    """Affiche les résultats des mesures"""
    st.subheader("📏 Résultats des Mesures")
    
    # Visualisation avec overlay
    try:
        engine = st.session_state.morpho_engine
        overlay_img = engine.generate_overlay_image(display_image, points, measurements)
        st.image(overlay_img, caption="Mesures annotées", use_container_width=True)
    except Exception as e:
        st.warning(f"Impossible de générer l'overlay: {str(e)}")
    
    # Tableau résultats
    col_res1, col_res2 = st.columns(2)
    
    with col_res1:
        st.markdown("**📐 Morphologie Générale**")
        metrics_general = {
            k: v for k, v in measurements.items() 
            if 'mamelle' not in k and 'score' not in k
        }
        for key, value in metrics_general.items():
            unit = " cm" if 'cm' in key else " kg" if 'kg' in key else ""
            st.metric(
                key.replace('_', ' ').title().replace('Cm', 'CM'), 
                f"{value}{unit}"
            )
    
    with col_res2:
        st.markdown("**🥛 Morphologie Mammelle**")
        if any('mamelle' in k for k in measurements.keys()):
            st.metric("Attachment", f"{measurements.get('attachment_mamelle_cm', 'N/A')} cm")
            st.metric("Profondeur", f"{measurements.get('profondeur_mamelle_cm', 'N/A')} cm")
            st.metric("Score Mamelle", f"{measurements.get('score_mamelle', 'N/A')}/10")
        else:
            st.info("Points mamelles non saisis")
        
        if 'poids_estime_kg' in measurements:
            st.metric("🎯 Poids Estimé", f"{measurements['poids_estime_kg']} kg")
    
    # Sauvegarde
    _render_save_section(measurements, db)

def _render_save_section(measurements: Dict[str, float], db: DatabaseManager):
    """Affiche la section de sauvegarde"""
    st.subheader("💾 Enregistrement")
    
    with st.form("save_morpho", clear_on_submit=True):
        animal_id = st.text_input("ID Animal", placeholder="FR123456789")
        notes = st.text_area("Notes", placeholder="Conditions de prise de vue, comportement...")
        
        col1, col2 = st.columns(2)
        with col1:
            submit = st.form_submit_button("💾 Sauvegarder en base")
        with col2:
            export = st.form_submit_button("📄 Exporter CSV")
        
        if submit and animal_id:
            try:
                db.execute(
                    """INSERT INTO mesures_morpho 
                       (brebis_id, date_mesure, longueur_cm, hauteur_garrot_cm, 
                        tour_poitrine_cm, circonf_canon_cm, score_mamelle, 
                        poids_estime_kg, notes, etalon_type, ppm_ratio)
                       VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                    (animal_id, date.today().isoformat(), 
                     measurements.get('longueur_corps_cm'),
                     measurements.get('hauteur_garrot_cm'),
                     measurements.get('tour_poitrine_cm'),
                     measurements.get('circonf_canon_cm'),
                     measurements.get('score_mamelle'),
                     measurements.get('poids_estime_kg'),
                     notes,
                     st.session_state.etalon_type.value,
                     st.session_state.ppm)
                )
                st.success(f"✅ Mesures enregistrées pour {animal_id}")
                
                # Réinitialisation partielle
                st.session_state.calibration_done = False
                st.session_state.last_measurements = None
                st.session_state.last_points = None
                
            except Exception as e:
                st.error(f"❌ Erreur sauvegarde: {str(e)}")
                logger.error(f"Erreur sauvegarde DB: {e}")
        elif submit:
            st.error("❌ ID animal requis")
        
        if export:
            # Export CSV
            df_export = pd.DataFrame([measurements])
            csv = df_export.to_csv(index=False)
            st.download_button(
                label="📥 Télécharger CSV",
                data=csv,
                file_name=f"mesures_morpho_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )

def _render_history_section(db: DatabaseManager):
    """Affiche l'historique des mesures"""
    st.subheader("📊 Historique des Mesures Scanner")
    
    try:
        df_morpho = db.fetch_df("""
            SELECT m.*, b.nom, b.race 
            FROM mesures_morpho m 
            LEFT JOIN brebis b ON m.brebis_id = b.identifiant_unique 
            ORDER BY m.date_mesure DESC 
            LIMIT 20
        """)
        
        if not df_morpho.empty:
            st.dataframe(df_morpho, use_container_width=True)
            
            # Graphique évolution si plusieurs mesures même animal
            animals_with_multi = df_morpho.groupby('brebis_id').size()
            multi_animals = animals_with_multi[animals_with_multi > 1].index.tolist()
            
            if multi_animals:
                selected_animal = st.selectbox("Voir évolution", multi_animals, key="select_animal_evo")
                df_evolution = df_morpho[df_morpho['brebis_id'] == selected_animal].sort_values('date_mesure')
                
                if len(df_evolution) > 1:
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=df_evolution['date_mesure'], 
                        y=df_evolution['poids_estime_kg'],
                        mode='lines+markers', 
                        name='Poids Estimé (kg)'
                    ))
                    fig.add_trace(go.Scatter(
                        x=df_evolution['date_mesure'], 
                        y=df_evolution['longueur_cm'],
                        mode='lines+markers', 
                        name='Longueur (cm)', 
                        yaxis='y2'
                    ))
                    
                    fig.update_layout(
                        title=f"Évolution morphométrique - {selected_animal}",
                        yaxis=dict(title="Poids (kg)"),
                        yaxis2=dict(
                            title="Longueur (cm)", 
                            overlaying='y', 
                            side='right'
                        ),
                        hovermode='x unified'
                    )
                    st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Aucune mesure enregistrée")
            
    except Exception as e:
        st.warning(f"Impossible de charger l'historique: {str(e)}")
        logger.error(f"Erreur chargement historique: {e}")

# ============================================================================
# DATABASE UPDATE (Ajout table mesures_morpho)
# ============================================================================

def init_morpho_schema(db: DatabaseManager):
    """Initialise le schéma de la base de données pour les mesures morpho"""
    schema = """
    CREATE TABLE IF NOT EXISTS brebis (
        identifiant_unique TEXT PRIMARY KEY,
        nom TEXT,
        race TEXT,
        date_naissance DATE,
        sexe TEXT,
        statut TEXT
    );
    
    CREATE TABLE IF NOT EXISTS mesures_morpho (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        brebis_id TEXT NOT NULL,
        date_mesure DATE NOT NULL,
        longueur_cm REAL,
        hauteur_garrot_cm REAL,
        tour_poitrine_cm REAL,
        circonf_canon_cm REAL,
        largeur_bassin_cm REAL,
        attachment_mamelle_cm REAL,
        profondeur_mamelle_cm REAL,
        score_mamelle REAL CHECK(score_mamelle BETWEEN 1 AND 10),
        poids_estime_kg REAL,
        indice_compacite REAL,
        notes TEXT,
        etalon_type TEXT,
        ppm_ratio REAL,  -- pixels per cm pour traçabilité
        image_ref BLOB,  -- Optionnel: stockage miniature
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY (brebis_id) REFERENCES brebis(identifiant_unique) ON DELETE CASCADE
    );
    
    CREATE INDEX IF NOT EXISTS idx_morpho_date ON mesures_morpho(date_mesure);
    CREATE INDEX IF NOT EXISTS idx_morpho_animal ON mesures_morpho(brebis_id);
    """
    
    try:
        db.execute("PRAGMA foreign_keys = ON")
        statements = [stmt.strip() for stmt in schema.split(';') if stmt.strip()]
        for statement in statements:
            try:
                db.execute(statement)
            except Exception as e:
                logger.warning(f"Exécution partielle du schéma: {statement[:50]}... - {e}")
        
        logger.info("Schéma morpho initialisé avec succès")
        return True
    except Exception as e:
        logger.error(f"Erreur initialisation schéma morpho: {e}")
        return False

# ============================================================================
# INTÉGRATION DANS LE ROUTING PRINCIPAL
# ============================================================================

def render_sidebar_updated() -> str:
    """Navigation avec module Scanner"""
    menu_items = {
        "📊 Dashboard": "Vue d'ensemble",
        "📝 Inscription": "Nouveaux animaux",
        "📷 Scanner IA": "Mesures 1m étalon",
        "🥛 Production": "Suivi laitier",
        "🩺 Santé": "Carnet sanitaire",
        "🧬 Génomique": "Analyse ADN",
        "🌾 Nutrition": "Rations",
        "⚙️ Admin": "Configuration"
    }
    
    return st.sidebar.radio(
        "Modules",
        list(menu_items.keys()),
        format_func=lambda x: f"{x} - {menu_items[x]}"
    )

# ============================================================================
# POINT D'ENTRÉE PRINCIPAL
# ============================================================================

def main():
    """Fonction principale de l'application"""
    st.set_page_config(
        page_title="EXPERT OVIN DZ PRO - Scanner IA",
        page_icon="🐑",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.sidebar.title("🐑 EXPERT OVIN DZ PRO")
    st.sidebar.markdown("---")
    
    # Navigation
    choice = render_sidebar_updated()
    
    # Initialisation DB
    try:
        with DatabaseManager() as db:
            # Initialisation du schéma (doit être fait une fois)
            init_morpho_schema(db)
            
            # Routage
            if choice == "📷 Scanner IA":
                render_scanner_ia(db)
            # ... autres choix de menu
            else:
                st.title(f"Module: {choice}")
                st.info("Module en développement")
                
    except Exception as e:
        st.error(f"❌ Erreur d'initialisation: {str(e)}")
        logger.critical(f"Erreur critique application: {e}", exc_info=True)

if __name__ == "__main__":
    main()
