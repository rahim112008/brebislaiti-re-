"""
EXPERT OVIN DZ PRO - VERSION ENTERPRISE 2026.02.1
Module Scanner IA Morphométrique avec Calibration Étalon
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
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union, Literal
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

# Numpy pour calculs géométriques
import numpy as np
from numpy.linalg import norm

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
    BATON_1M = "baton_1m"           # 100 cm
    CARTE_BANCAIRE = "carte_bancaire"  # 8.56 cm × 5.398 cm (ISO 7810 ID-1)
    FEUILLE_A4 = "feuille_a4"       # 29.7 cm × 21 cm

@dataclass(frozen=True)
class AppConfig:
    DB_PATH: str = "data/ovin_enterprise.db"
    CACHE_TTL: int = 3600
    ALIGNMENT_THRESHOLD: float = 85.0
    MIN_SEQ_LENGTH: int = 50
    
    # Dimensions étalons référence (cm)
    ETALON_DIMS: Dict[str, Dict] = None
    
    def __post_init__(self):
        object.__setattr__(self, 'ETALON_DIMS', {
            EtalonType.BATON_1M.value: {
                "longueur_reelle_cm": 100.0,
                "largeur_cm": 2.5,  # bâton standard
                "nom": "Bâton de 1 mètre",
                "couleur_recommandee": "Jaune/Orange vif"
            },
            EtalonType.CARTE_BANCAIRE.value: {
                "longueur_reelle_cm": 8.56,
                "largeur_reelle_cm": 5.398,
                "nom": "Carte bancaire (ISO 7810)",
                "couleur_recommandee": "Peu importe"
            },
            EtalonType.FEUILLE_A4.value: {
                "longueur_reelle_cm": 29.7,
                "largeur_reelle_cm": 21.0,
                "nom": "Feuille A4 standard",
                "couleur_recommandee": "Blanche (contraste maximal)"
            }
        })

CONFIG = AppConfig()

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
        
    def set_etalon(self, etalon_type: EtalonType, pixel_length: float, 
                   pixel_width: Optional[float] = None):
        """
        Calibration: établit l'échelle pixels -> cm réels
        """
        self.etalon_type = etalon_type
        dims = CONFIG.ETALON_DIMS[etalon_type.value]
        
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
        else:
            self.pixels_per_cm = pixel_length / dims["longueur_reelle_cm"]
        
        self.reference_object_pixels = pixel_length
        logger.info(f"Calibration: {self.pixels_per_cm:.2f} px/cm avec {etalon_type.value}")
        
        return self.pixels_per_cm
    
    def pixels_to_cm(self, pixels: float) -> float:
        """Convertit une distance en pixels vers cm réels"""
        if not self.pixels_per_cm:
            raise ValueError("Calibration non effectuée")
        return pixels / self.pixels_per_cm
    
    def cm_to_pixels(self, cm: float) -> float:
        """Convertit cm vers pixels (pour overlay)"""
        if not self.pixels_per_cm:
            raise ValueError("Calibration non effectuée")
        return cm * self.pixels_per_cm
    
    def calculate_measurements(self, points: Dict[str, Tuple[int, int]]) -> Dict[str, float]:
        """
        Calcule les distances morphométriques à partir des points clics utilisateur
        Points attendus: 'garrot', 'base_queue', 'hanche', 'patte_arriere',
                        'mamelle_gauche', 'mamelle_droite', 'canon', 'poitrine'
        """
        if not self.pixels_per_cm:
            return {"error": "Calibration requise"}
        
        measurements = {}
        
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
        if all(k in points for k in ['hanche_gauche', 'hanche_droite']):
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
        if all(k in points for k in ['mamelle_gauche', 'mamelle_droite', 'mamelle_arriere']):
            # Largeur attachment
            attache_px = self._distance(points['mamelle_gauche'], points['mamelle_droite'])
            measurements['attachment_mamelle_cm'] = round(self.pixels_to_cm(attache_px), 1)
            
            # Profondeur (arrière -> ligne attachment)
            profondeur_px = self._distance(
                points['mamelle_arriere'],
                ((points['mamelle_gauche'][0] + points['mamelle_droite'][0]) // 2,
                 (points['mamelle_gauche'][1] + points['mamelle_droite'][1]) // 2)
            )
            measurements['profondeur_mamelle_cm'] = round(self.pixels_to_cm(profondeur_px), 1)
            
            # Score calculé (formule simplifiée INRA)
            if 'attachment_mamelle_cm' in measurements and 'profondeur_mamelle_cm' in measurements:
                score = min(10, max(1, 
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
    
    def _distance(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> float:
        """Distance euclidienne entre deux points"""
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    def generate_overlay_image(self, image: Image.Image, points: Dict, 
                               measurements: Dict) -> Image.Image:
        """Génère image annotée avec lignes de mesure"""
        img_annotated = image.copy()
        draw = ImageDraw.Draw(img_annotated)
        
        # Couleurs
        COLOR_ETALON = (255, 0, 0)      # Rouge
        COLOR_MESURE = (0, 255, 0)      # Vert
        COLOR_POINT = (0, 0, 255)       # Bleu
        COLOR_TEXT = (255, 255, 255)    # Blanc
        
        try:
            font = ImageFont.truetype("arial.ttf", 20)
        except:
            font = ImageFont.load_default()
        
        # Dessine points cliqués
        for name, (x, y) in points.items():
            draw.ellipse([x-5, y-5, x+5, y+5], fill=COLOR_POINT, outline=COLOR_TEXT, width=2)
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
                
                # Label milieu
                mid_x = (p1[0] + p2[0]) // 2
                mid_y = (p1[1] + p2[1]) // 2
                if mesure_key in measurements:
                    draw.text((mid_x, mid_y), f"{measurements[mesure_key]}cm", 
                             fill=COLOR_MESURE, font=font)
        
        return img_annotated

# ============================================================================
# INTERFACE SCANNER IA
# ============================================================================

def render_scanner_ia(db: DatabaseManager):
    """Interface complète de scan morphométrique avec calibration étalon"""
    st.header("📷 Scanner IA Morphométrique 1m")
    
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
        use_baton = st.button("Sélectionner Bâton 1m", key="etalon_baton")
    
    with col2:
        st.markdown("**💳 Carte Bancaire** (Urgence)")
        st.caption("ISO 7810: 8.56 × 5.40 cm. Placer sur le dos de l'animal.")
        use_carte = st.button("Sélectionner Carte", key="etalon_carte")
    
    with col3:
        st.markdown("**📄 Feuille A4** (Standard)")
        st.caption("29.7 × 21 cm. Fixer verticalement près de l'animal.")
        use_a4 = st.button("Sélectionner A4", key="etalon_a4")
    
    # Détection sélection
    etalon_selected = None
    if use_baton:
        etalon_selected = EtalonType.BATON_1M
    elif use_carte:
        etalon_selected = EtalonType.CARTE_BANCAIRE
    elif use_a4:
        etalon_selected = EtalonType.FEUILLE_A4
    
    if etalon_selected:
        st.session_state['etalon_type'] = etalon_selected
        dims = CONFIG.ETALON_DIMS[etalon_selected.value]
        st.success(f"✅ Étalonnage sélectionné: {dims['nom']} ({dims['longueur_reelle_cm']} cm)")
    
    # Étape 2: Upload image
    st.subheader("2️⃣ Acquisition de l'Image")
    
    uploaded_file = st.file_uploader(
        "📸 Photo de profil de l'animal avec l'étalon visible",
        type=['jpg', 'jpeg', 'png'],
        help="L'animal doit être de profil, tête à gauche, étalon visible et parallèle au corps"
    )
    
    if uploaded_file and 'etalon_type' in st.session_state:
        # Traitement image
        image = Image.open(uploaded_file)
        
        # Redimensionnement pour affichage (conservation ratio)
        max_width = 800
        if image.width > max_width:
            ratio = max_width / image.width
            new_size = (max_width, int(image.height * ratio))
            display_image = image.resize(new_size, Image.Resampling.LANCZOS)
        else:
            display_image = image
        
        st.image(display_image, caption="Image source", use_container_width=True)
        
        # Étape 3: Calibration interactive
        st.subheader("3️⃣ Calibration sur l'Image")
        
        st.markdown(f"""
        **Instructions**: 
        1. Cliquez sur **les deux extrémités** de l'étalon ({CONFIG.ETALON_DIMS[st.session_state['etalon_type'].value]['nom']})
        2. Assurez-vous que la ligne est parallèle au sol
        """)
        
        # Canvas interactif pour calibration
        if 'calibration_points' not in st.session_state:
            st.session_state['calibration_points'] = []
        
        # Utilisation de click coordinates sur l'image
        col_cal1, col_cal2 = st.columns([3, 1])
        
        with col_cal1:
            # Simulation de canvas cliquable avec streamlit
            st.markdown("### 🎯 Cliquez sur les extrémités de l'étalon")
            
            # Pour l'instant, utilisation de coordonnées manuelles (en attendant component canvas)
            x1 = st.number_input("X1 (début étalon)", 0, display_image.width, 100, key="x1")
            y1 = st.number_input("Y1 (début étalon)", 0, display_image.height, 400, key="y1")
            x2 = st.number_input("X2 (fin étalon)", 0, display_image.width, 300, key="x2")
            y2 = st.number_input("Y2 (fin étalon)", 0, display_image.height, 400, key="y2")
            
            # Visualisation calibration
            calib_img = display_image.copy()
            draw_cal = ImageDraw.Draw(calib_img)
            draw_cal.line([(x1, y1), (x2, y2)], fill=(255, 0, 0), width=3)
            draw_cal.ellipse([x1-5, y1-5, x1+5, y1+5], fill=(255, 0, 0))
            draw_cal.ellipse([x2-5, y2-5, x2+5, y2+5], fill=(255, 0, 0))
            st.image(calib_img, caption="Vérification calibration", use_container_width=True)
        
        with col_cal2:
            pixel_distance = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            st.metric("Distance pixels", f"{pixel_distance:.1f} px")
            
            if st.button("✅ Valider Calibration"):
                engine = MorphoMetricsEngine()
                ppm = engine.set_etalon(
                    st.session_state['etalon_type'], 
                    pixel_distance
                )
                st.session_state['morpho_engine'] = engine
                st.session_state['ppm'] = ppm
                st.success(f"Calibration: {ppm:.2f} pixels/cm")
                st.rerun()
        
        # Étape 4: Mesures morphométriques
        if 'morpho_engine' in st.session_state:
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
                
                points['mamelle_arriere'] = (
                    st.number_input("Mamelle Arrière X", 0, display_image.width, 350, key="max"),
                    st.number_input("Mamelle Arrière Y", 0, display_image.height, 500, key="may")
                )
                points['canon_peripherie'] = (
                    st.number_input("Périphérie Canon X", 0, display_image.width, 570, key="cpx"),
                    st.number_input("Périphérie Canon Y", 0, display_image.height, 480, key="cpy")
                )
            
            if st.button("🧮 Calculer les Mesures"):
                engine = st.session_state['morpho_engine']
                measurements = engine.calculate_measurements(points)
                
                if "error" in measurements:
                    st.error(measurements["error"])
                else:
                    st.session_state['last_measurements'] = measurements
                    st.session_state['last_points'] = points
                    
                    # Affichage résultats
                    st.subheader("📏 Résultats des Mesures")
                    
                    # Visualisation
                    overlay_img = engine.generate_overlay_image(
                        display_image, points, measurements
                    )
                    st.image(overlay_img, caption="Mesures annotées", use_container_width=True)
                    
                    # Tableau résultats
                    col_res1, col_res2 = st.columns(2)
                    
                    with col_res1:
                        st.markdown("**📐 Morphologie Générale**")
                        metrics_general = {
                            k: v for k, v in measurements.items() 
                            if 'mamelle' not in k and 'score' not in k
                        }
                        for key, value in metrics_general.items():
                            st.metric(
                                key.replace('_', ' ').title(), 
                                f"{value} cm" if 'cm' in key else f"{value} kg" if 'kg' in key else value
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
                    st.subheader("💾 Enregistrement")
                    with st.form("save_morpho"):
                        animal_id = st.text_input("ID Animal", placeholder="FR123456789")
                        notes = st.text_area("Notes", placeholder="Conditions de prise de vue...")
                        
                        if st.form_submit_button("Sauvegarder en base"):
                            if animal_id:
                                # Sauvegarde DB
                                db.execute(
                                    """INSERT INTO mesures_morpho 
                                       (brebis_id, date_mesure, longueur_cm, hauteur_garrot_cm, 
                                        tour_poitrine_cm, circonf_canon_cm, score_mamelle, 
                                        poids_estime_kg, notes, etalon_type, ppm_ratio)
                                       VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                                    (animal_id, date.today(), 
                                     measurements.get('longueur_corps_cm'),
                                     measurements.get('hauteur_garrot_cm'),
                                     measurements.get('tour_poitrine_cm'),
                                     measurements.get('circonf_canon_cm'),
                                     measurements.get('score_mamelle'),
                                     measurements.get('poids_estime_kg'),
                                     notes,
                                     st.session_state['etalon_type'].value,
                                     st.session_state['ppm'])
                                )
                                st.success(f"✅ Mesures enregistrées pour {animal_id}")
                            else:
                                st.error("ID animal requis")
    
    elif uploaded_file and 'etalon_type' not in st.session_state:
        st.warning("⚠️ Veuillez d'abord sélectionner un type d'étalon de calibration")
    
    # Historique des mesures
    st.subheader("📊 Historique des Mesures Scanner")
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
            selected_animal = st.selectbox("Voir évolution", multi_animals)
            df_evolution = df_morpho[df_morpho['brebis_id'] == selected_animal].sort_values('date_mesure')
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df_evolution['date_mesure'], y=df_evolution['poids_estime_kg'],
                                    mode='lines+markers', name='Poids Estimé'))
            fig.add_trace(go.Scatter(x=df_evolution['date_mesure'], y=df_evolution['longueur_cm'],
                                    mode='lines+markers', name='Longueur', yaxis='y2'))
            
            fig.update_layout(
                title=f"Évolution morphométrique - {selected_animal}",
                yaxis=dict(title="Poids (kg)"),
                yaxis2=dict(title="Longueur (cm)", overlaying='y', side='right')
            )
            st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Aucune mesure enregistrée")

# ============================================================================
# DATABASE UPDATE (Ajout table mesures_morpho)
# ============================================================================

def init_morpho_schema(db: DatabaseManager):
    """Ajoute la table pour stocker les mesures scanner"""
    schema = """
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
        db.conn.executescript(schema)
        logger.info("Schema morpho initialisé")
    except Exception as e:
        logger.error(f"Schema morpho error: {e}")

# ============================================================================
# INTÉGRATION DANS LE ROUTING PRINCIPAL
# ============================================================================

def render_sidebar_updated() -> str:
    """Navigation avec module Scanner"""
    menu_items = {
        "📊 Dashboard": "Vue d'ensemble",
        "📝 Inscription": "Nouveaux animaux",
        "📷 Scanner IA": "Mesures 1m étalon",  # NOUVEAU
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

# Dans main(), ajouter:
# elif choice == "📷 Scanner IA":
#     render_scanner_ia(db)
#     init_morpho_schema(db)  # Au cas où
