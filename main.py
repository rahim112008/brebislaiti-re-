"""
EXPERT OVIN DZ PRO - VERSION MASTER 2026.02.09 AVEC MESURES MAMMAIRES
Système Intégral : Phénotypage, Bio-Informatique (GWAS Pro, PLINK), 
Accouplement IA & Alignement Génomique de Référence (Ensembl/EBI).
Reconnaissance d'Animaux par IA avec Mesures Morphométriques + Mamelles.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import re
import requests
import cv2
import mediapipe as mp
from datetime import datetime, date, timedelta
import random
import tempfile
from PIL import Image
import json
import base64
import io
import math
from scipy.spatial import distance
from shapely.geometry import Polygon, Point

# ============================================================================
# 1. DATABASE MASTER
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_pro.db"):
        self.db_path = db_path
        if not os.path.exists('data'): os.makedirs('data')
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row

    def execute_query(self, query: str, params: tuple = ()):
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            if "duplicate column name" not in str(e).lower():
                st.error(f"Erreur SQL: {e}")
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, sexe TEXT, age_type TEXT, age_valeur REAL,
            hauteur REAL, longueur REAL, tour_poitrine REAL, 
            largeur_bassin REAL, long_bassin REAL, circ_canon REAL,
            note_mamelle INTEGER, attaches_mamelle TEXT, poids REAL, created_at DATE,
            image_path TEXT, race_detectee TEXT, confiance_race REAL,
            -- Nouvelles mesures mammaires
            largeur_mamelle REAL, hauteur_mamelle REAL, profondeur_mamelle REAL,
            distance_tetines REAL, diametre_tetine REAL, symetrie_mamelle REAL,
            volume_mamelle_estime REAL, score_morpho_mamelle INTEGER
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_controle DATE,
            quantite_lait REAL, tb REAL, tp REAL, cellules INTEGER,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS genomique (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, marqueur TEXT, zygotie TEXT, impact TEXT, date_test DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_eponge DATE, date_mise_bas_prevue DATE, statut TEXT
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )""",
        """CREATE TABLE IF NOT EXISTS images_animals (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id TEXT, image_path TEXT, date_capture DATE,
            race_detectee TEXT, confiance REAL, mesures TEXT,
            mesures_mamelles TEXT, -- Nouveau: mesures mammaires en JSON
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS references_mesures (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            session_id TEXT, type_reference TEXT, 
            longueur_reelle_cm REAL, longueur_pixels REAL,
            date_utilisation DATE, angle_vue TEXT
        )"""
    ]
    for table_sql in tables: 
        db.execute_query(table_sql)
    
    # Ajouter les colonnes si elles n'existent pas
    columns_to_add = [
        ("largeur_mamelle", "REAL"), ("hauteur_mamelle", "REAL"),
        ("profondeur_mamelle", "REAL"), ("distance_tetines", "REAL"),
        ("diametre_tetine", "REAL"), ("symetrie_mamelle", "REAL"),
        ("volume_mamelle_estime", "REAL"), ("score_morpho_mamelle", "INTEGER")
    ]
    
    for col_name, col_type in columns_to_add:
        try:
            db.execute_query(f"ALTER TABLE brebis ADD COLUMN {col_name} {col_type}")
        except:
            pass  # Colonne déjà existante

def seed_data_demo(db: DatabaseManager):
    races = ["Ouled Djellal", "Rembi", "Hamra", "Lacaune"]
    markers = ["CAST", "DGAT1", "PrP", "GDF8"]
    for i in range(1, 16):
        uid = f"DZ-2026-{100+i}"
        sexe = "Mâle" if i > 12 else "Femelle"
        race = random.choice(races)
        db.execute_query("""INSERT OR IGNORE INTO brebis 
            (identifiant_unique, nom, race, sexe, age_type, age_valeur, hauteur, longueur, tour_poitrine, circ_canon, note_mamelle, poids, created_at) 
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (uid, f"Animal_{i}", race, sexe, "Années", random.randint(2, 5), 75, 80, 95, 8.5, random.randint(4, 9), random.randint(55, 85), date.today()))

        for m in markers:
            db.execute_query("INSERT OR IGNORE INTO genomique (brebis_id, marqueur, zygotie, impact, date_test) VALUES (?,?,?,?,?)",
                             (uid, m, random.choice(["Homozygote", "Hétérozygote", "Absent"]), "Auto-Généré", date.today()))

# ============================================================================
# 2. SYSTÈME DE CALIBRATION AVEC RÉFÉRENCES MULTIPLES
# ============================================================================

class CalibrationSystem:
    """Système de calibration avec objets de référence standards"""
    
    REFERENCE_OBJECTS = {
        "baton_1m": {
            "longueur_cm": 100.0,
            "description": "Bâton ou mètre de 1 mètre",
            "couleur_recommandee": "rouge/jaune pour détection facile",
            "tolérance_cm": 0.5
        },
        "feuille_a4": {
            "longueur_cm": 29.7,  # Hauteur A4
            "largeur_cm": 21.0,   # Largeur A4
            "description": "Feuille A4 standard",
            "ratio": 1.414,  # √2
            "tolérance_cm": 0.2
        },
        "carte_bancaire": {
            "longueur_cm": 8.56,
            "largeur_cm": 5.398,
            "description": "Carte bancaire (format CR80)",
            "ratio": 1.586,
            "tolérance_cm": 0.05
        },
        "piece_10da": {
            "diametre_cm": 2.7,
            "description": "Pièce de 10 DA algérienne",
            "tolérance_cm": 0.05
        }
    }
    
    @staticmethod
    def detect_reference_object(image_path, object_type):
        """
        Détecte un objet de référence dans l'image
        """
        image = cv2.imread(image_path)
        if image is None:
            return None
        
        h, w = image.shape[:2]
        
        # Convertir en HSV pour détection de couleur
        hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
        
        # Définir les plages de couleurs pour différents objets
        color_ranges = {
            "baton_1m": [
                (np.array([0, 100, 100]), np.array([10, 255, 255])),  # Rouge
                (np.array([20, 100, 100]), np.array([30, 255, 255]))  # Jaune
            ],
            "feuille_a4": [
                (np.array([0, 0, 200]), np.array([180, 30, 255]))  # Blanc
            ],
            "carte_bancaire": [
                (np.array([100, 150, 0]), np.array([140, 255, 255]))  # Bleu/Vert
            ]
        }
        
        if object_type not in color_ranges:
            # Détection par contours pour tout objet
            gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            blurred = cv2.GaussianBlur(gray, (5, 5), 0)
            edges = cv2.Canny(blurred, 50, 150)
            
            contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            if contours:
                # Filtrer les contours par taille et forme
                for contour in contours:
                    area = cv2.contourArea(contour)
                    if 1000 < area < 50000:  # Taille raisonnable pour référence
                        x, y, w_cont, h_cont = cv2.boundingRect(contour)
                        
                        # Calculer le ratio hauteur/largeur
                        ratio = h_cont / w_cont if w_cont > 0 else 0
                        
                        # Vérifier si c'est probablement un objet de référence
                        ref_info = CalibrationSystem.REFERENCE_OBJECTS.get(object_type, {})
                        expected_ratio = ref_info.get('ratio', None)
                        
                        if expected_ratio and 0.8*expected_ratio < ratio < 1.2*expected_ratio:
                            return {
                                "bbox": (x, y, w_cont, h_cont),
                                "center": (x + w_cont/2, y + h_cont/2),
                                "pixel_length": max(w_cont, h_cont),
                                "confidence": 0.7
                            }
            
            return None
        
        # Détection par couleur pour objets spécifiques
        masks = []
        for lower, upper in color_ranges[object_type]:
            mask = cv2.inRange(hsv, lower, upper)
            masks.append(mask)
        
        if masks:
            combined_mask = cv2.bitwise_or(masks[0], masks[1]) if len(masks) > 1 else masks[0]
            
            # Appliquer une morphologie pour nettoyer
            kernel = np.ones((5,5), np.uint8)
            combined_mask = cv2.morphologyEx(combined_mask, cv2.MORPH_CLOSE, kernel)
            combined_mask = cv2.morphologyEx(combined_mask, cv2.MORPH_OPEN, kernel)
            
            contours, _ = cv2.findContours(combined_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            if contours:
                # Prendre le plus grand contour
                largest_contour = max(contours, key=cv2.contourArea)
                x, y, w_cont, h_cont = cv2.boundingRect(largest_contour)
                
                return {
                    "bbox": (x, y, w_cont, h_cont),
                    "center": (x + w_cont/2, y + h_cont/2),
                    "pixel_length": max(w_cont, h_cont),
                    "confidence": 0.8
                }
        
        return None
    
    @staticmethod
    def calculate_scale_factor(object_type, pixel_length, reference_info=None):
        """
        Calcule l'échelle pixels → cm
        """
        ref_data = CalibrationSystem.REFERENCE_OBJECTS.get(object_type, {})
        
        if object_type == "baton_1m":
            real_length_cm = ref_data["longueur_cm"]
        elif object_type == "feuille_a4":
            # Utiliser la hauteur (29.7 cm) par défaut
            real_length_cm = ref_data["longueur_cm"]
        elif object_type == "carte_bancaire":
            real_length_cm = ref_data["longueur_cm"]
        elif object_type == "piece_10da":
            real_length_cm = ref_data["diametre_cm"]
        else:
            # Si info personnalisée fournie
            real_length_cm = reference_info.get("longueur_cm", 100.0)
        
        if pixel_length > 0:
            scale_cm_per_pixel = real_length_cm / pixel_length
            return scale_cm_per_pixel
        
        return None
    
    @staticmethod
    def draw_reference_markers(image_path, detected_objects):
        """
        Dessine des marqueurs sur les objets de référence détectés
        """
        image = cv2.imread(image_path)
        if image is None:
            return None
        
        for obj_type, detection in detected_objects.items():
            if detection:
                x, y, w, h = detection["bbox"]
                color = (0, 255, 0)  # Vert pour référence
                
                # Dessiner le rectangle
                cv2.rectangle(image, (x, y), (x+w, y+h), color, 2)
                
                # Ajouter le texte
                label = f"{obj_type}: {detection['pixel_length']}px"
                cv2.putText(image, label, (x, y-10), 
                           cv2.FONT_HERSHEY_SIMPLEX, 0.5, color, 2)
                
                # Dessiner le centre
                cx, cy = detection["center"]
                cv2.circle(image, (int(cx), int(cy)), 5, (0, 0, 255), -1)
        
        return image

# ============================================================================
# 3. ANALYSE MAMMAIRE POUR OVINS
# ============================================================================

class UdderAnalyzer:
    """Analyseur spécifique pour la mamelle des ovins"""
    
    @staticmethod
    def detect_udder_region(image_path, scale_cm_per_pixel=None):
        """
        Détecte et analyse la région mammaire
        """
        try:
            image = cv2.imread(image_path)
            if image is None:
                return None
            
            h, w = image.shape[:2]
            
            # 1. Détection de la zone ventrale (bas de l'abdomen)
            # Convertir en HSV pour détection de peau/poils
            hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
            
            # Plage pour poils ovins (brun, blanc, noir)
            lower_skin = np.array([0, 10, 60])
            upper_skin = np.array([20, 150, 255])
            
            mask = cv2.inRange(hsv, lower_skin, upper_skin)
            
            # 2. Trouver les contours de la zone ventrale
            contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            if not contours:
                # Fallback: utiliser la partie inférieure centrale de l'image
                udder_zone = {
                    "bbox": (w//4, 3*h//4, w//2, h//4),
                    "center": (w//2, 7*h//8),
                    "confidence": 0.3,
                    "teats_detected": False
                }
                return udder_zone
            
            # Filtrer les contours pour trouver la zone la plus basse (mamelle)
            lowest_contour = None
            lowest_y = 0
            
            for contour in contours:
                _, y_cont, _, h_cont = cv2.boundingRect(contour)
                bottom_y = y_cont + h_cont
                
                # Vérifier si c'est dans le tiers inférieur de l'image
                if bottom_y > 2*h/3 and bottom_y > lowest_y:
                    lowest_contour = contour
                    lowest_y = bottom_y
            
            if lowest_contour is not None:
                x, y, w_cont, h_cont = cv2.boundingRect(lowest_contour)
                
                # Ajuster pour zone mammaire (plus large que haute)
                udder_height = max(30, min(h_cont, h//10))  # Limiter la hauteur
                udder_width = min(w_cont * 2, w//3)  # Élargir un peu
                
                udder_zone = {
                    "bbox": (x, y, udder_width, udder_height),
                    "center": (x + udder_width//2, y + udder_height//2),
                    "confidence": 0.6,
                    "teats_detected": False,
                    "contour_area": cv2.contourArea(lowest_contour)
                }
                
                # 3. Détection des tétines (si visible)
                udder_zone = UdderAnalyzer._detect_teats(image, udder_zone, scale_cm_per_pixel)
                
                return udder_zone
            
            return None
            
        except Exception as e:
            st.warning(f"Erreur détection mamelle: {e}")
            return None
    
    @staticmethod
    def _detect_teats(image, udder_zone, scale_cm_per_pixel):
        """
        Détecte les tétines dans la zone mammaire
        """
        x, y, w_zone, h_zone = udder_zone["bbox"]
        
        # Extraire la ROI (Region of Interest)
        roi = image[y:y+h_zone, x:x+w_zone]
        
        if roi.size == 0:
            return udder_zone
        
        # Convertir en niveaux de gris
        gray_roi = cv2.cvtColor(roi, cv2.COLOR_BGR2GRAY)
        
        # Appliquer un flou pour réduire le bruit
        blurred = cv2.GaussianBlur(gray_roi, (5, 5), 0)
        
        # Détection de cercles (tétines)
        circles = cv2.HoughCircles(
            blurred, 
            cv2.HOUGH_GRADIENT, 
            dp=1, 
            minDist=20,
            param1=50, 
            param2=30,
            minRadius=5,
            maxRadius=30
        )
        
        teats = []
        if circles is not None:
            circles = np.uint16(np.around(circles))
            
            for circle in circles[0, :]:
                # Convertir les coordonnées relatives en absolues
                cx_abs = x + circle[0]
                cy_abs = y + circle[1]
                radius = circle[2]
                
                teats.append({
                    "center": (cx_abs, cy_abs),
                    "radius": radius,
                    "diameter_px": radius * 2
                })
        
        udder_zone["teats"] = teats
        udder_zone["teats_detected"] = len(teats) >= 2
        
        # Si échelle disponible, convertir en cm
        if scale_cm_per_pixel and udder_zone["teats_detected"]:
            for teat in teats:
                teat["diameter_cm"] = teat["diameter_px"] * scale_cm_per_pixel
        
        return udder_zone
    
    @staticmethod
    def calculate_udder_measurements(udder_zone, scale_cm_per_pixel=None):
        """
        Calcule les mesures mammaires détaillées
        """
        if not udder_zone:
            return None
        
        x, y, w_zone, h_zone = udder_zone["bbox"]
        
        mesures = {
            "largeur_mamelle_px": w_zone,
            "hauteur_mamelle_px": h_zone,
            "profondeur_estimee_px": min(w_zone, h_zone) * 0.6,  # Estimation
            "symetrie_score": 0.5,
            "volume_estime_ml": 0,
            "nombre_tetines": len(udder_zone.get("teats", [])),
            "teats": udder_zone.get("teats", [])
        }
        
        # Calculer la distance entre tétines si au moins 2 détectées
        if mesures["nombre_tetines"] >= 2:
            teats = udder_zone["teats"]
            # Distance entre les deux tétines les plus centrales
            if len(teats) >= 2:
                # Trier par position X
                teats_sorted = sorted(teats, key=lambda t: t["center"][0])
                # Prendre les deux plus proches du centre
                center_x = x + w_zone/2
                teats_centered = sorted(teats, 
                                      key=lambda t: abs(t["center"][0] - center_x))
                
                if len(teats_centered) >= 2:
                    t1, t2 = teats_centered[:2]
                    dist_px = math.sqrt(
                        (t2["center"][0] - t1["center"][0])**2 + 
                        (t2["center"][1] - t1["center"][1])**2
                    )
                    mesures["distance_tetines_px"] = dist_px
                    
                    # Score de symétrie (plus proche de 1 = plus symétrique)
                    if t1["radius"] > 0 and t2["radius"] > 0:
                        mesures["symetrie_score"] = min(t1["radius"], t2["radius"]) / max(t1["radius"], t2["radius"])
        
        # Convertir en cm si échelle disponible
        if scale_cm_per_pixel:
            mesures["largeur_mamelle_cm"] = w_zone * scale_cm_per_pixel
            mesures["hauteur_mamelle_cm"] = h_zone * scale_cm_per_pixel
            mesures["profondeur_mamelle_cm"] = mesures["profondeur_estimee_px"] * scale_cm_per_pixel
            
            if "distance_tetines_px" in mesures:
                mesures["distance_tetines_cm"] = mesures["distance_tetines_px"] * scale_cm_per_pixel
            
            # Estimer le volume mammaire (formule simplifiée)
            # Volume ≈ largeur × hauteur × profondeur × facteur
            volume_cm3 = (mesures["largeur_mamelle_cm"] * 
                         mesures["hauteur_mamelle_cm"] * 
                         mesures["profondeur_mamelle_cm"] * 0.5)
            mesures["volume_estime_ml"] = round(volume_cm3, 1)
            
            # Calculer le diamètre moyen des tétines
            if mesures["teats"]:
                total_diameter = sum(t.get("diameter_cm", t["diameter_px"] * scale_cm_per_pixel) 
                                   for t in mesures["teats"])
                mesures["diametre_tetine_moyen_cm"] = total_diameter / len(mesures["teats"])
        
        # Calculer le score morphologique (1-9)
        mesures["score_morphologique"] = UdderAnalyzer._calculate_udder_score(mesures)
        
        return mesures
    
    @staticmethod
    def _calculate_udder_score(mesures):
        """
        Calcule un score morphologique de mamelle (1-9)
        Basé sur les standards ovins
        """
        score = 5  # Moyen par défaut
        
        if "largeur_mamelle_cm" in mesures:
            largeur = mesures["largeur_mamelle_cm"]
            hauteur = mesures.get("hauteur_mamelle_cm", largeur * 0.7)
            
            # Score basé sur ratio largeur/hauteur (idéal ~1.5)
            ratio = largeur / hauteur if hauteur > 0 else 1
            if 1.3 <= ratio <= 1.7:
                score += 2
            elif 1.1 <= ratio <= 1.9:
                score += 1
            
            # Score basé sur symétrie
            symetrie = mesures.get("symetrie_score", 0.5)
            if symetrie > 0.9:
                score += 2
            elif symetrie > 0.8:
                score += 1
            
            # Score basé sur nombre de tétines (idéal: 2)
            nb_tetines = mesures.get("nombre_tetines", 0)
            if nb_tetines == 2:
                score += 1
            elif nb_tetines > 2:
                score -= 1  # Tétines supplémentaires moins désirables
        
        # Limiter entre 1 et 9
        return max(1, min(9, int(round(score))))
    
    @staticmethod
    def draw_udder_analysis(image_path, udder_zone, measurements):
        """
        Dessine l'analyse mammaire sur l'image
        """
        image = cv2.imread(image_path)
        if image is None:
            return None
        
        if udder_zone:
            # Dessiner la zone mammaire
            x, y, w, h = udder_zone["bbox"]
            cv2.rectangle(image, (x, y), (x+w, y+h), (0, 255, 255), 2)  # Jaune
            
            # Dessiner les tétines
            if "teats" in udder_zone:
                for teat in udder_zone["teats"]:
                    cx, cy = teat["center"]
                    radius = teat["radius"]
                    cv2.circle(image, (cx, cy), radius, (0, 0, 255), 2)  # Rouge
                    cv2.circle(image, (cx, cy), 3, (255, 0, 0), -1)  # Centre bleu
            
            # Ajouter les mesures
            if measurements:
                y_text = y - 10 if y > 30 else y + h + 20
                text = f"Mamelle: {measurements.get('largeur_mamelle_cm', 'N/A')}cm"
                cv2.putText(image, text, (x, y_text), 
                           cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 255), 2)
        
        return image

# ============================================================================
# 4. MOTEURS IA AMÉLIORÉS
# ============================================================================

class AnimalRecognitionAPI:
    """API de reconnaissance d'animaux et mesures morphométriques"""
    
    @staticmethod
    def detect_animal_species(image_file, api_key=None):
        """
        Détecte l'espèce animale dans l'image
        """
        # (Le code reste identique à la version précédente)
        try:
            if api_key:
                return AnimalRecognitionAPI._google_vision_detect(image_file, api_key)
            else:
                return AnimalRecognitionAPI._local_animal_detect(image_file)
        except Exception as e:
            st.error(f"Erreur reconnaissance: {e}")
            return {"species": "Ovin", "confidence": 0.7, "breeds": ["Ouled Djellal"]}

    # ... (autres méthodes inchangées)

class MorphometricAnalyzer:
    """Analyse morphométrique améliorée avec calibration"""
    
    def __init__(self):
        self.mp_pose = mp.solutions.pose
        self.pose = self.mp_pose.Pose(
            static_image_mode=True,
            model_complexity=1,
            enable_segmentation=False,
            min_detection_confidence=0.5
        )
        self.calibration_system = CalibrationSystem()
    
    def analyze_image(self, image_path, reference_scale=None, detected_references=None):
        """
        Analyse améliorée avec calibration par référence
        """
        try:
            image = cv2.imread(image_path)
            if image is None:
                raise ValueError("Impossible de lire l'image")
            
            h, w, _ = image.shape
            
            # Si des références sont détectées, calculer l'échelle précise
            scale_cm_per_pixel = reference_scale
            if detected_references:
                scale_cm_per_pixel = self._calculate_precise_scale(detected_references)
            
            # Détection de la pose
            image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
            results = self.pose.process(image_rgb)
            
            mesures = {
                "hauteur_pixels": 0,
                "longueur_pixels": 0,
                "tour_poitrine_pixels": 0,
                "largeur_bassin_pixels": 0,
                "longueur_bassin_pixels": 0,
                "circ_canon_pixels": 0,
                "confidence": 0,
                "landmarks_detected": 0
            }
            
            if results.pose_landmarks:
                landmarks = results.pose_landmarks.landmark
                mesures["landmarks_detected"] = len(landmarks)
                
                # Calcul des mesures avec points anatomiques précis
                mesures.update(self._calculate_precise_measurements(landmarks, w, h))
                mesures["confidence"] = min(0.85, len(landmarks) / 15)
            
            # Alternative: détection par contours
            if mesures["landmarks_detected"] < 5:
                mesures.update(self._contour_analysis(image))
            
            # Convertir en cm si échelle disponible
            if scale_cm_per_pixel:
                mesures.update(self._convert_to_cm(mesures, scale_cm_per_pixel))
            
            return mesures
            
        except Exception as e:
            st.warning(f"Analyse morphométrique limitée: {e}")
            return self._get_default_measurements()
    
    def _calculate_precise_scale(self, detected_references):
        """
        Calcule l'échelle moyenne à partir de plusieurs références
        """
        scales = []
        
        for ref_type, detection in detected_references.items():
            if detection:
                ref_data = CalibrationSystem.REFERENCE_OBJECTS.get(ref_type, {})
                if ref_data:
                    real_length = ref_data.get("longueur_cm") or ref_data.get("diametre_cm", 100)
                    pixel_length = detection["pixel_length"]
                    
                    if pixel_length > 0:
                        scale = real_length / pixel_length
                        scales.append(scale)
        
        if scales:
            # Prendre la médiane pour éviter les outliers
            return np.median(scales)
        
        return None
    
    def _calculate_precise_measurements(self, landmarks, image_width, image_height):
        """
        Calcule des mesures précises à partir des landmarks
        """
        mesures = {}
        
        # Points clés (indices MediaPipe adaptés pour animaux)
        points = {}
        for idx, lm in enumerate(landmarks):
            points[idx] = (lm.x * image_width, lm.y * image_height)
        
        # 1. Hauteur au garrot (withers height)
        # Utiliser le point le plus haut du dos (landmark 8 ou 7)
        if 8 in points and 27 in points:  # Garrot à sol via patte
            y_garrot = points[8][1]
            y_sol = points[27][1]
            mesures["hauteur_pixels"] = abs(y_garrot - y_sol)
        
        # 2. Longueur du corps (épaule à fesse)
        if 11 in points and 23 in points:  # Épaule droite à hanche droite
            x_epaule = points[11][0]
            x_hanche = points[23][0]
            mesures["longueur_pixels"] = abs(x_hanche - x_epaule)
        
        # 3. Largeur de bassin (distance entre hanches)
        if 23 in points and 24 in points:
            largeur = abs(points[23][0] - points[24][0])
            mesures["largeur_bassin_pixels"] = largeur
        
        # 4. Longueur de bassin (hanche à fesse)
        if 23 in points and 25 in points:
            longueur = abs(points[23][0] - points[25][0])
            mesures["longueur_bassin_pixels"] = longueur
        
        # 5. Tour de poitrine (estimation géométrique)
        if all(k in points for k in [7, 11, 12]):
            # Largeur entre épaules
            largeur_ep = abs(points[11][0] - points[12][0])
            # Profondeur poitrine (sternum à dos)
            profondeur = abs(points[7][1] - points[8][1])
            # Estimation circonférence: 2×(largeur+profondeur)×0.7
            mesures["tour_poitrine_pixels"] = 2 * (largeur_ep + profondeur) * 0.7
        
        # 6. Circonférence du canon (estimation)
        if 27 in points and 28 in points:
            # Distance entre les canons × π
            distance_canons = abs(points[27][0] - points[28][0])
            mesures["circ_canon_pixels"] = distance_canons * math.pi * 0.3
        
        return mesures
    
    def _convert_to_cm(self, mesures, scale_cm_per_pixel):
        """Convertit toutes les mesures de pixels en cm"""
        cm_measures = {}
        
        for key, value in mesures.items():
            if key.endswith("_pixels") and isinstance(value, (int, float)):
                cm_key = key.replace("_pixels", "_cm")
                cm_measures[cm_key] = value * scale_cm_per_pixel
        
        return cm_measures
    
    # ... (autres méthodes inchangées)

# ============================================================================
# 5. MODULE RECONNAISSANCE ANIMALE AMÉLIORÉ
# ============================================================================

def animal_recognition_module(db):
    """Module amélioré avec calibration et analyse mammaire"""
    st.title("📷 Reconnaissance IA des Animaux & Mesures Morphométriques")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("1. Capture d'Image avec Références")
        
        # Option 1: Téléchargement d'image
        uploaded_file = st.file_uploader("📤 Télécharger une photo d'animal", 
                                        type=['jpg', 'jpeg', 'png', 'bmp'])
        
        # Option 2: Webcam
        use_camera = st.checkbox("Utiliser la caméra")
        camera_image = None
        if use_camera:
            camera_image = st.camera_input("Prendre une photo")
        
        image_to_process = uploaded_file or camera_image
        
        if image_to_process:
            # Sauvegarder l'image temporairement
            with tempfile.NamedTemporaryFile(delete=False, suffix='.jpg') as tmp_file:
                if uploaded_file:
                    tmp_file.write(uploaded_file.getvalue())
                else:
                    tmp_file.write(camera_image.getvalue())
                tmp_path = tmp_file.name
            
            # Afficher l'image
            st.image(image_to_process, caption="Image à analyser", use_column_width=True)
            
            # Section CALIBRATION avec références
            st.subheader("2. 🔧 Calibration avec Objets de Référence")
            
            calibration_cols = st.columns(4)
            
            with calibration_cols[0]:
                use_baton = st.checkbox("📏 Bâton 1m", value=True)
            with calibration_cols[1]:
                use_a4 = st.checkbox("📄 Feuille A4")
            with calibration_cols[2]:
                use_carte = st.checkbox("💳 Carte bancaire")
            with calibration_cols[3]:
                custom_ref = st.checkbox("📐 Référence perso.")
            
            custom_length = None
            if custom_ref:
                custom_length = st.number_input("Longueur référence (cm)", 
                                               min_value=1.0, max_value=500.0, 
                                               value=100.0, step=0.1)
            
            # Information sur l'animal
            st.subheader("3. 🐑 Informations de l'Animal")
            animal_id = st.text_input("Identifiant de l'animal", 
                                     value=f"DZ-{datetime.now().year}-{random.randint(1000, 9999)}")
            nom_animal = st.text_input("Nom (optionnel)", value="")
            sexe_animal = st.selectbox("Sexe", ["Femelle", "Mâle", "Inconnu"])
            
            # Options d'analyse
            st.subheader("4. ⚙️ Options d'Analyse")
            col_opt1, col_opt2 = st.columns(2)
            with col_opt1:
                analyze_udder = st.checkbox("🔍 Analyser la mamelle", value=True)
            with col_opt2:
                estimate_weight = st.checkbox("⚖️ Estimer le poids", value=True)
            
            if st.button("🚀 Lancer l'Analyse Complète", type="primary"):
                with st.spinner("Analyse en cours..."):
                    # Initialiser les analyseurs
                    recognizer = AnimalRecognitionAPI()
                    morpho_analyzer = MorphometricAnalyzer()
                    udder_analyzer = UdderAnalyzer()
                    calibration_system = CalibrationSystem()
                    
                    # Étape 1: Détection des objets de référence
                    detected_references = {}
                    if use_baton or use_a4 or use_carte:
                        st.info("📏 Calibration avec références...")
                        
                        if use_baton:
                            detected_references["baton_1m"] = calibration_system.detect_reference_object(
                                tmp_path, "baton_1m")
                        
                        if use_a4:
                            detected_references["feuille_a4"] = calibration_system.detect_reference_object(
                                tmp_path, "feuille_a4")
                        
                        if use_carte:
                            detected_references["carte_bancaire"] = calibration_system.detect_reference_object(
                                tmp_path, "carte_bancaire")
                    
                    # Étape 2: Calcul de l'échelle
                    scale_cm_per_pixel = None
                    if detected_references:
                        # Afficher les références détectées
                        ref_image = calibration_system.draw_reference_markers(tmp_path, detected_references)
                        if ref_image is not None:
                            st.image(ref_image, caption="Objets de référence détectés", use_column_width=True)
                        
                        # Calculer l'échelle
                        for ref_type, detection in detected_references.items():
                            if detection:
                                scale = calibration_system.calculate_scale_factor(ref_type, 
                                                                                 detection["pixel_length"])
                                if scale:
                                    scale_cm_per_pixel = scale
                                    st.success(f"Échelle calculée: 1px = {scale:.4f} cm")
                                    break
                    
                    # Étape 3: Reconnaissance de l'espèce
                    st.info("🔍 Reconnaissance de l'espèce...")
                    detection_results = recognizer.detect_animal_species(tmp_path)
                    
                    # Étape 4: Analyse morphométrique
                    st.info("📐 Calcul des mesures corporelles...")
                    morpho_results = morpho_analyzer.analyze_image(tmp_path, scale_cm_per_pixel, detected_references)
                    
                    # Étape 5: Analyse mammaire (si femelle et option activée)
                    udder_results = None
                    udder_image = None
                    
                    if analyze_udder and sexe_animal == "Femelle":
                        st.info("🥛 Analyse de la mamelle...")
                        
                        # Détecter la mamelle
                        udder_zone = udder_analyzer.detect_udder_region(tmp_path, scale_cm_per_pixel)
                        
                        if udder_zone:
                            # Calculer les mesures mammaires
                            udder_results = udder_analyzer.calculate_udder_measurements(
                                udder_zone, scale_cm_per_pixel)
                            
                            # Dessiner l'analyse
                            udder_image = udder_analyzer.draw_udder_analysis(
                                tmp_path, udder_zone, udder_results)
                            
                            if udder_image is not None:
                                st.image(udder_image, caption="Analyse de la mamelle", use_column_width=True)
                    
                    # Étape 6: Estimation du poids
                    weight_estimate = None
                    if estimate_weight and "tour_poitrine_cm" in morpho_results:
                        # Formule pour ovins: Poids (kg) = (Tour de poitrine² × Longueur) / 10800
                        tour = morpho_results.get("tour_poitrine_cm", 90)
                        longueur = morpho_results.get("longueur_cm", 80)
                        weight_estimate = (tour**2 * longueur) / 10800
                    
                    # Étape 7: Affichage des résultats
                    st.success("✅ Analyse terminée !")
                    
                    # Afficher les résultats dans la colonne de droite
                    with col2:
                        st.subheader("📋 Résultats")
                        
                        # Informations générales
                        st.metric("Animal ID", animal_id)
                        st.metric("Sexe", sexe_animal)
                        
                        # Espèce détectée
                        if detection_results:
                            st.metric("Espèce", detection_results.get("species", "Ovin"))
                            st.metric("Confiance", f"{detection_results.get('confidence', 0)*100:.1f}%")
                        
                        # Mesures corporelles principales
                        st.subheader("📏 Mesures Corporelles")
                        if "hauteur_cm" in morpho_results:
                            st.metric("Hauteur", f"{morpho_results['hauteur_cm']:.1f} cm")
                        if "longueur_cm" in morpho_results:
                            st.metric("Longueur", f"{morpho_results['longueur_cm']:.1f} cm")
                        if "tour_poitrine_cm" in morpho_results:
                            st.metric("Tour poitrine", f"{morpho_results['tour_poitrine_cm']:.1f} cm")
                        
                        # Mesures mammaires si disponibles
                        if udder_results:
                            st.subheader("🥛 Mesures Mammaires")
                            if "largeur_mamelle_cm" in udder_results:
                                st.metric("Largeur mamelle", f"{udder_results['largeur_mamelle_cm']:.1f} cm")
                            if "hauteur_mamelle_cm" in udder_results:
                                st.metric("Hauteur mamelle", f"{udder_results['hauteur_mamelle_cm']:.1f} cm")
                            if "distance_tetines_cm" in udder_results:
                                st.metric("Distance tétines", f"{udder_results['distance_tetines_cm']:.1f} cm")
                            if "score_morphologique" in udder_results:
                                score = udder_results["score_morphologique"]
                                color = "🟢" if score >= 7 else "🟡" if score >= 5 else "🔴"
                                st.metric("Score mamelle", f"{color} {score}/9")
                        
                        # Poids estimé
                        if weight_estimate:
                            st.subheader("⚖️ Estimation Poids")
                            st.metric("Poids estimé", f"{weight_estimate:.1f} kg")
                    
                    # Étape 8: Enregistrement dans la base de données
                    st.subheader("💾 Sauvegarde des Résultats")
                    
                    # Préparer les données JSON
                    mesures_json = json.dumps({
                        "morphometrie": morpho_results,
                        "detection": detection_results,
                        "udder": udder_results,
                        "calibration": {
                            "scale_cm_per_pixel": scale_cm_per_pixel,
                            "references_detected": bool(detected_references)
                        },
                        "date_analyse": datetime.now().isoformat()
                    })
                    
                    try:
                        # Sauvegarder dans la table images_animals
                        db.execute_query(
                            """INSERT INTO images_animals 
                            (brebis_id, image_path, date_capture, race_detectee, confiance, mesures, mesures_mamelles) 
                            VALUES (?, ?, ?, ?, ?, ?, ?)""",
                            (animal_id, tmp_path, date.today(), 
                             detection_results.get("primary_breed", "Inconnue"),
                             detection_results.get("confidence", 0.5),
                             mesures_json,
                             json.dumps(udder_results) if udder_results else None)
                        )
                        
                        # Préparer les données pour la table brebis
                        race_connue = detection_results.get("primary_breed") or "Ouled Djellal"
                        
                        # Données de base
                        brebis_data = {
                            'identifiant_unique': animal_id,
                            'nom': nom_animal or f"Animal_{animal_id}",
                            'race': race_connue,
                            'race_detectee': race_connue,
                            'confiance_race': detection_results.get("confidence", 0.7),
                            'sexe': sexe_animal,
                            'image_path': tmp_path,
                            'created_at': date.today()
                        }
                        
                        # Ajouter les mesures corporelles
                        if morpho_results:
                            for key in ['hauteur_cm', 'longueur_cm', 'tour_poitrine_cm', 
                                       'largeur_bassin_cm', 'longueur_bassin_cm', 'circ_canon_cm']:
                                if key in morpho_results:
                                    db_key = key.replace('_cm', '')
                                    brebis_data[db_key] = morpho_results[key]
                        
                        # Ajouter les mesures mammaires
                        if udder_results:
                            udder_mapping = {
                                'largeur_mamelle_cm': 'largeur_mamelle',
                                'hauteur_mamelle_cm': 'hauteur_mamelle',
                                'profondeur_mamelle_cm': 'profondeur_mamelle',
                                'distance_tetines_cm': 'distance_tetines',
                                'diametre_tetine_moyen_cm': 'diametre_tetine',
                                'symetrie_score': 'symetrie_mamelle',
                                'volume_estime_ml': 'volume_mamelle_estime',
                                'score_morphologique': 'score_morpho_mamelle'
                            }
                            
                            for src_key, dst_key in udder_mapping.items():
                                if src_key in udder_results:
                                    brebis_data[dst_key] = udder_results[src_key]
                        
                        # Ajouter le poids estimé
                        if weight_estimate:
                            brebis_data['poids'] = weight_estimate
                        
                        # Vérifier si l'animal existe déjà
                        existing = db.fetch_all_as_df(
                            "SELECT * FROM brebis WHERE identifiant_unique = ?",
                            (animal_id,)
                        )
                        
                        if existing.empty:
                            # Créer une nouvelle entrée
                            columns = ', '.join(brebis_data.keys())
                            placeholders = ', '.join(['?' for _ in brebis_data])
                            values = tuple(brebis_data.values())
                            
                            query = f"INSERT INTO brebis ({columns}) VALUES ({placeholders})"
                            db.execute_query(query, values)
                            st.success(f"✅ Nouvel animal {animal_id} enregistré !")
                        else:
                            # Mettre à jour l'entrée existante
                            set_clause = ', '.join([f"{k} = ?" for k in brebis_data.keys()])
                            values = tuple(brebis_data.values()) + (animal_id,)
                            
                            query = f"UPDATE brebis SET {set_clause} WHERE identifiant_unique = ?"
                            db.execute_query(query, values)
                            st.success(f"✅ Animal {animal_id} mis à jour !")
                        
                        # Afficher un récapitulatif
                        st.balloons()
                        
                        # Récapitulatif détaillé
                        with st.expander("📊 Récapitulatif détaillé"):
                            st.json(brebis_data)
                        
                    except Exception as e:
                        st.error(f"Erreur lors de l'enregistrement: {e}")
    
    # Section historique
    st.divider()
    st.subheader("📊 Historique des Analyses")
    
    if st.button("📈 Afficher l'historique complet"):
        historique = db.fetch_all_as_df(
            """SELECT i.brebis_id, i.race_detectee, i.confiance, i.date_capture,
                      b.hauteur, b.longueur, b.tour_poitrine, b.score_morpho_mamelle
               FROM images_animals i
               LEFT JOIN brebis b ON i.brebis_id = b.identifiant_unique
               ORDER BY i.date_capture DESC LIMIT 20"""
        )
        
        if not historique.empty:
            st.dataframe(historique)
            
            # Graphiques
            col_chart1, col_chart2 = st.columns(2)
            
            with col_chart1:
                fig1 = px.pie(historique, names='race_detectee', 
                             title='Distribution des races détectées')
                st.plotly_chart(fig1)
            
            with col_chart2:
                if 'score_morpho_mamelle' in historique.columns:
                    fig2 = px.histogram(historique, x='score_morpho_mamelle',
                                       title='Distribution des scores mammaires',
                                       nbins=9, range_x=[1, 9])
                    st.plotly_chart(fig2)
            
            # Corrélations
            st.subheader("🔗 Corrélations entre mesures")
            corr_data = historique[['hauteur', 'longueur', 'tour_poitrine']].dropna()
            if not corr_data.empty:
                fig3 = px.scatter_matrix(corr_data, title="Matrice de corrélation")
                st.plotly_chart(fig3)
        else:
            st.info("Aucune analyse enregistrée.")

# ============================================================================
# 6. INTERFACE UTILISATEUR PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🧬")
    
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
        st.session_state.morpho_analyzer = MorphometricAnalyzer()
        st.session_state.udder_analyzer = UdderAnalyzer()
    
    db = st.session_state.db
    morpho_analyzer = st.session_state.morpho_analyzer
    udder_analyzer = st.session_state.udder_analyzer
    ia = AIEngine()
    api_web = WebBioAPI()

    st.sidebar.title("🐑 Bio-Master v2026")
    if st.sidebar.button("🚀 Charger Données Démo Pro"):
        seed_data_demo(db)
        st.sidebar.success("Base initialisée !")
    
    # Instructions de calibration
    with st.sidebar.expander("📏 Guide de Calibration"):
        st.info("""
        **Pour des mesures précises:**
        
        1. **Bâton 1m**: Peint en rouge/jaune
        2. **Feuille A4**: Blanche, placée verticalement
        3. **Carte bancaire**: Standard (8.56×5.398 cm)
        
        Placez la référence à la même distance que l'animal.
        """)
    
    # Mettre à jour le menu
    menu = [
        "📊 Dashboard Élite", "🧬 Génomique & Alignement", "📈 GWAS & PLINK Pro", 
        "🌐 Recherche Bio-Web", "⚤ Accouplement IA", "🥛 Contrôle Laitier", 
        "📷 Reconnaissance Animale IA", "🤰 Gestation IA", "🌾 Nutrition Solo", 
        "📈 Statistiques Mammaires"  # Nouveau menu
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- MODULE : RECONNAISSANCE ANIMALE IA ---
    if choice == "📷 Reconnaissance Animale IA":
        animal_recognition_module(db)

    # --- NOUVEAU MODULE : STATISTIQUES MAMMAIRES ---
    elif choice == "📈 Statistiques Mammaires":
        st.title("📈 Statistiques Mammaires Avancées")
        
        # Récupérer les données mammaires
        query = """
        SELECT identifiant_unique, race, sexe,
               largeur_mamelle, hauteur_mamelle, distance_tetines,
               score_morpho_mamelle, volume_mamelle_estime
        FROM brebis 
        WHERE largeur_mamelle IS NOT NULL 
          AND sexe = 'Femelle'
        ORDER BY score_morpho_mamelle DESC
        """
        
        df_mamelles = db.fetch_all_as_df(query)
        
        if not df_mamelles.empty:
            st.success(f"✅ {len(df_mamelles)} animaux avec données mammaires")
            
            # Tableau des meilleures mamelles
            st.subheader("🏆 Top 10 des meilleures mamelles")
            top_10 = df_mamelles.head(10)
            st.dataframe(top_10)
            
            # Graphiques
            col1, col2 = st.columns(2)
            
            with col1:
                # Distribution des scores
                fig1 = px.histogram(df_mamelles, x='score_morpho_mamelle',
                                   title='Distribution des scores mammaires (1-9)',
                                   nbins=9, color='race')
                st.plotly_chart(fig1)
            
            with col2:
                # Relation largeur-hauteur
                fig2 = px.scatter(df_mamelles, x='largeur_mamelle', y='hauteur_mamelle',
                                 color='score_morpho_mamelle', size='volume_mamelle_estime',
                                 title='Relation Largeur-Hauteur de mamelle',
                                 hover_data=['identifiant_unique', 'race'])
                st.plotly_chart(fig2)
            
            # Corrélations par race
            st.subheader("📊 Corrélations par race")
            races = df_mamelles['race'].unique()
            
            for race in races:
                df_race = df_mamelles[df_mamelles['race'] == race]
                if len(df_race) > 3:
                    st.write(f"**{race}** (n={len(df_race)})")
                    
                    # Calculer les moyennes
                    avg_width = df_race['largeur_mamelle'].mean()
                    avg_height = df_race['hauteur_mamelle'].mean()
                    avg_score = df_race['score_morpho_mamelle'].mean()
                    
                    col_avg1, col_avg2, col_avg3 = st.columns(3)
                    col_avg1.metric("Largeur moyenne", f"{avg_width:.1f} cm")
                    col_avg2.metric("Hauteur moyenne", f"{avg_height:.1f} cm")
                    col_avg3.metric("Score moyen", f"{avg_score:.1f}/9")
        
        else:
            st.info("Aucune donnée mammaire disponible. Utilisez le module de reconnaissance pour analyser des femelles.")

    # --- AUTRES MODULES (inchangés) ---
    elif choice == "🧬 Génomique & Alignement":
        # ... (code existant)
        pass
    elif choice == "📊 Dashboard Élite":
        # ... (code existant)
        pass
    # ... (autres modules)

# ============================================================================
# INSTALLATION DES DEPENDANCES
# ============================================================================

def install_requirements():
    """Fonction pour installer les dépendances nécessaires"""
    requirements = """
    # Pour la reconnaissance animale et mesures morphométriques:
    opencv-python>=4.8.0
    mediapipe>=0.10.0
    Pillow>=10.0.0
    numpy>=1.24.0
    scipy>=1.11.0
    shapely>=2.0.0
    
    # Pour l'API Google Vision (optionnel):
    google-cloud-vision>=3.0.0
    
    # Dépendances existantes:
    streamlit>=1.28.0
    pandas>=2.0.0
    plotly>=5.17.0
    requests>=2.31.0
    sqlite3
    """
    return requirements

if __name__ == "__main__":
    # Afficher les instructions d'installation
    with st.sidebar.expander("📦 Installation"):
        st.code(install_requirements(), language="bash")
        st.info("""
        **Pour installer:**
        ```bash
        pip install opencv-python mediapipe Pillow scipy shapely
        ```
        """)
    
    main()
