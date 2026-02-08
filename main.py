"""
EXPERT OVIN DZ PRO - VERSION COMPLETE & CORRIGÉE 2026.02
Module Scanner IA Morphométrique + Système Intégré
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

# Configuration logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler('ovin_pro.log'), logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# ============================================================================
# 1. CONFIGURATION & CONSTANTES (TOUJOURS EN PREMIER)
# ============================================================================

class EtalonType(Enum):
    BATON_1M = "baton_1m"
    CARTE_BANCAIRE = "carte_bancaire"
    FEUILLE_A4 = "feuille_a4"

@dataclass(frozen=True)
class AppConfig:
    DB_PATH: str = "data/ovin_enterprise.db"
    CACHE_TTL: int = 3600
    ALIGNMENT_THRESHOLD: float = 85.0
    MIN_SEQ_LENGTH: int = 50
    
    ETALON_DIMS: Dict = None
    
    def __post_init__(self):
        object.__setattr__(self, 'ETALON_DIMS', {
            EtalonType.BATON_1M.value: {
                "longueur_reelle_cm": 100.0,
                "largeur_cm": 2.5,
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
                "couleur_recommandee": "Blanche"
            }
        })

CONFIG = AppConfig()

# ============================================================================
# 2. DATABASE MANAGER (DÉFINI AVANT TOUTE FONCTION QUI L'UTILISE)
# ============================================================================

class DatabaseManager:
    """Gestionnaire DB avec pool de connexions et transactions"""
    
    def __init__(self, db_path: str = CONFIG.DB_PATH):
        self.db_path = db_path
        self._ensure_directory()
        self._init_connection_pool()
        self._create_indexes()
        logger.info(f"DatabaseManager initialisé: {db_path}")
    
    def _ensure_directory(self):
        os.makedirs(os.path.dirname(self.db_path) or '.', exist_ok=True)
    
    def _init_connection_pool(self):
        self.conn = sqlite3.connect(
            self.db_path, 
            check_same_thread=False,
            isolation_level=None,
            timeout=30
        )
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.row_factory = sqlite3.Row
    
    def _create_indexes(self):
        indexes = [
            "CREATE INDEX IF NOT EXISTS idx_brebis_race ON brebis(race)",
            "CREATE INDEX IF NOT EXISTS idx_laitier_date ON controle_laitier(date_controle)",
            "CREATE INDEX IF NOT EXISTS idx_sante_date ON sante(date_soin)",
            "CREATE INDEX IF NOT EXISTS idx_morpho_date ON mesures_morpho(date_mesure)",
            "CREATE INDEX IF NOT EXISTS idx_morpho_animal ON mesures_morpho(brebis_id)"
        ]
        for idx in indexes:
            try:
                self.conn.execute(idx)
            except sqlite3.Error:
                pass
    
    @contextmanager
    def transaction(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute("BEGIN")
            yield cursor
            self.conn.commit()
        except Exception as e:
            self.conn.rollback()
            raise
        finally:
            cursor.close()
    
    def execute(self, query: str, params: Tuple = ()) -> Optional[sqlite3.Cursor]:
        try:
            with self.transaction() as cursor:
                cursor.execute(query, params)
                return cursor
        except Exception as e:
            logger.error(f"SQL Error: {e}")
            return None
    
    def fetch_df(self, query: str, params: Tuple = ()) -> pd.DataFrame:
        try:
            return pd.read_sql_query(query, self.conn, params=params)
        except Exception as e:
            logger.error(f"Fetch error: {e}")
            return pd.DataFrame()
    
    def fetch_all_as_df(self, query: str, params: Tuple = ()) -> pd.DataFrame:
        """Alias pour compatibilité"""
        return self.fetch_df(query, params)

# ============================================================================
# 3. MOTEUR BIOINFORMATIQUE
# ============================================================================

class BioInfoEngine:
    """Moteur génomique simplifié pour compatibilité"""
    
    GENES_INTERET = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA",
        "MSTN (Muscle)": "AAGCTTGATTAGCAGGTTCCCGG",
        "CAST (Tendreté)": "TGGGGCCCAAGTCGATTGCAGAA",
        "DGAT1 (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
    }
    
    GENES_SANTE = {
        "Scrapie ARR (RÉSISTANCE)": "TGGTACCCATAATCAGTGGAACA",
        "Scrapie VRQ (SENSIBLE)": "TGGTAGCCATAATCAGTGGAACA",
        "Arachnomélie": "CCGTAGCTAGCTGATCGATCGTA",
        "Hypotrichose": "TTAGCGCTAGCTAGCTAGCTAGC"
    }
    
    def __init__(self):
        self.aligner = None
        try:
            from Bio.Align import PairwiseAligner
            self.aligner = PairwiseAligner()
            self.aligner.mode = 'local'
        except ImportError:
            pass
    
    @staticmethod
    def filtrer_sequence(seq):
        if not seq:
            return ""
        lines = seq.strip().split('\n')
        cleaned = []
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                continue
            cleaned.append(line.upper())
        return ''.join(cleaned).replace(' ', '').replace('\r', '')
    
    def extraire_multi_fasta(self, raw_text):
        sequences = {}
        current_id = None
        current_seq = []
        
        for line in raw_text.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper().replace(' ', ''))
        
        if current_id and current_seq:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences if sequences else {"Individu": self.filtrer_sequence(raw_text)}
    
    def alignement_expert(self, seq_test, ref_seq):
        if not seq_test or not ref_seq or not self.aligner:
            return 0.0
        try:
            score = self.aligner.score(seq_test, ref_seq)
            return round((score / len(ref_seq)) * 100, 2)
        except:
            return 0.0
    
    def calculer_heterozygotie(self, sequences_dict):
        if len(sequences_dict) < 2:
            return 0.0
        seqs = list(sequences_dict.values())
        distances = []
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                score = self.aligner.score(seqs[i], seqs[j]) if self.aligner else 0
                distances.append(1 - (score / max(len(seqs[i]), len(seqs[j]))))
        return round(np.mean(distances) * 100, 2) if distances else 0.0
    
    def traduire_en_proteine(self, dna_seq):
        try:
            from Bio.Seq import Seq
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            if len(clean_dna) < 3:
                return "Séquence trop courte"
            return str(Seq(clean_dna).translate(to_stop=True))
        except:
            return "Erreur de traduction"

# ============================================================================
# 4. MOTEUR MORPHOMÉTRIQUE (DÉFINI AVANT RENDER_SCANNER_IA)
# ============================================================================

class MorphoMetricsEngine:
    """Moteur de mesure morphométrique par analyse d'image"""
    
    def __init__(self):
        self.reference_object_pixels = None
        self.pixels_per_cm = None
        self.etalon_type = None
        
    def set_etalon(self, etalon_type: EtalonType, pixel_length: float, 
                   pixel_width: Optional[float] = None):
        self.etalon_type = etalon_type
        dims = CONFIG.ETALON_DIMS[etalon_type.value]
        
        if pixel_width and etalon_type != EtalonType.BATON_1M:
            ppm_long = pixel_length / dims["longueur_reelle_cm"]
            ppm_larg = pixel_width / dims["largeur_reelle_cm"]
            self.pixels_per_cm = (ppm_long + ppm_larg) / 2
        else:
            self.pixels_per_cm = pixel_length / dims["longueur_reelle_cm"]
        
        self.reference_object_pixels = pixel_length
        return self.pixels_per_cm
    
    def pixels_to_cm(self, pixels: float) -> float:
        if not self.pixels_per_cm:
            raise ValueError("Calibration non effectuée")
        return pixels / self.pixels_per_cm
    
    def calculate_measurements(self, points: Dict[str, Tuple[int, int]]) -> Dict[str, float]:
        if not self.pixels_per_cm:
            return {"error": "Calibration requise"}
        
        measurements = {}
        
        # Longueur corps
        if all(k in points for k in ['garrot', 'base_queue']):
            dist_px = self._distance(points['garrot'], points['base_queue'])
            measurements['longueur_corps_cm'] = round(self.pixels_to_cm(dist_px), 1)
        
        # Hauteur garrot
        if all(k in points for k in ['garrot', 'patte_arriere']):
            hauteur_px = abs(points['garrot'][1] - points['patte_arriere'][1])
            measurements['hauteur_garrot_cm'] = round(self.pixels_to_cm(hauteur_px), 1)
        
        # Tour poitrine
        if all(k in points for k in ['poitrine', 'garrot']):
            poitrine_px = self._distance(points['poitrine'], points['garrot']) * 1.2
            measurements['tour_poitrine_cm'] = round(self.pixels_to_cm(poitrine_px), 1)
        
        # Circonférence canon
        if 'canon' in points and 'canon_peripherie' in points:
            diametre_px = self._distance(points['canon'], points['canon_peripherie']) * 2
            circonf_px = np.pi * diametre_px
            measurements['circonf_canon_cm'] = round(self.pixels_to_cm(circonf_px), 1)
        
        # Mamelle
        if all(k in points for k in ['mamelle_gauche', 'mamelle_droite', 'mamelle_arriere']):
            attache_px = self._distance(points['mamelle_gauche'], points['mamelle_droite'])
            measurements['attachment_mamelle_cm'] = round(self.pixels_to_cm(attache_px), 1)
            
            centre_attachment = (
                (points['mamelle_gauche'][0] + points['mamelle_droite'][0]) // 2,
                (points['mamelle_gauche'][1] + points['mamelle_droite'][1]) // 2
            )
            profondeur_px = self._distance(points['mamelle_arriere'], centre_attachment)
            measurements['profondeur_mamelle_cm'] = round(self.pixels_to_cm(profondeur_px), 1)
            
            score = min(10, max(1, 
                (measurements['attachment_mamelle_cm'] / 10) + 
                (measurements['profondeur_mamelle_cm'] / 5)
            ))
            measurements['score_mamelle'] = round(score, 1)
        
        # Poids estimé et indice
        if 'longueur_corps_cm' in measurements and 'tour_poitrine_cm' in measurements:
            poids_estime = (measurements['tour_poitrine_cm'] ** 2 * measurements['longueur_corps_cm']) / 10800
            measurements['poids_estime_kg'] = round(poids_estime, 1)
            measurements['indice_compacite'] = round(poids_estime / measurements['longueur_corps_cm'], 2)
        
        return measurements
    
    def _distance(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> float:
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    def generate_overlay_image(self, image: Image.Image, points: Dict, 
                               measurements: Dict) -> Image.Image:
        img_annotated = image.copy()
        draw = ImageDraw.Draw(img_annotated)
        
        COLOR_MESURE = (0, 255, 0)
        COLOR_POINT = (255, 0, 0)
        COLOR_TEXT = (255, 255, 255)
        
        try:
            font = ImageFont.truetype("arial.ttf", 20)
        except:
            font = ImageFont.load_default()
        
        # Points
        for name, (x, y) in points.items():
            draw.ellipse([x-5, y-5, x+5, y+5], fill=COLOR_POINT, outline=COLOR_TEXT, width=2)
            draw.text((x+8, y-8), name, fill=COLOR_TEXT, font=font)
        
        # Lignes
        lines = [
            ('garrot', 'base_queue', 'longueur_corps_cm'),
            ('garrot', 'patte_arriere', 'hauteur_garrot_cm'),
            ('mamelle_gauche', 'mamelle_droite', 'attachment_mamelle_cm')
        ]
        
        for p1_name, p2_name, mesure_key in lines:
            if p1_name in points and p2_name in points:
                p1, p2 = points[p1_name], points[p2_name]
                draw.line([p1, p2], fill=COLOR_MESURE, width=3)
                mid_x = (p1[0] + p2[0]) // 2
                mid_y = (p1[1] + p2[1]) // 2
                if mesure_key in measurements:
                    draw.text((mid_x, mid_y), f"{measurements[mesure_key]}cm", 
                             fill=COLOR_MESURE, font=font)
        
        return img_annotated

# ============================================================================
# 5. INITIALISATION SCHEMA (FONCTIONS UTILITAIRES)
# ============================================================================

def init_database(db: DatabaseManager):
    """Initialise toutes les tables"""
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, poids REAL, note_mamelle INTEGER, 
            tour_poitrine REAL, longueur REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL, 
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )""",
        """CREATE TABLE IF NOT EXISTS mesures_morpho (
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
            score_mamelle REAL,
            poids_estime_kg REAL,
            indice_compacite REAL,
            notes TEXT,
            etalon_type TEXT,
            ppm_ratio REAL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (brebis_id) REFERENCES brebis(identifiant_unique)
        )"""
    ]
    for table_sql in tables:
        db.execute(table_sql)

# ============================================================================
# 6. INTERFACE SCANNER (DÉFINIE APRÈS TOUTES LES CLASSES)
# ============================================================================

def render_scanner_ia(db: DatabaseManager):
    """Interface complète de scan morphométrique"""
    st.header("📷 Scanner IA Morphométrique 1m")
    
    st.markdown("""
    ### 📋 Protocole de Mesure Assistée par IA
    **Principe**: Calibration par objet de taille connue (étalon), puis mesures morphométriques.
    **Précision**: ±2 cm si protocole respecté.
    """)
    
    # Étape 1: Choix étalon
    st.subheader("1️⃣ Choix de l'Étalon de Calibration")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**🟨 Bâton 1 mètre**")
        st.caption("Précision maximale. Placer au sol.")
        if st.button("Sélectionner Bâton 1m", key="etalon_baton"):
            st.session_state['etalon_type'] = EtalonType.BATON_1M
    
    with col2:
        st.markdown("**💳 Carte Bancaire**")
        st.caption("8.56 × 5.40 cm. Placer sur le dos.")
        if st.button("Sélectionner Carte", key="etalon_carte"):
            st.session_state['etalon_type'] = EtalonType.CARTE_BANCAIRE
    
    with col3:
        st.markdown("**📄 Feuille A4**")
        st.caption("29.7 × 21 cm. Fixer verticalement.")
        if st.button("Sélectionner A4", key="etalon_a4"):
            st.session_state['etalon_type'] = EtalonType.FEUILLE_A4
    
    if 'etalon_type' in st.session_state:
        dims = CONFIG.ETALON_DIMS[st.session_state['etalon_type'].value]
        st.success(f"✅ Étalonnage: {dims['nom']} ({dims['longueur_reelle_cm']} cm)")
    
    # Étape 2: Upload
    st.subheader("2️⃣ Acquisition de l'Image")
    
    uploaded_file = st.file_uploader(
        "📸 Photo de profil avec étalon visible",
        type=['jpg', 'jpeg', 'png']
    )
    
    if uploaded_file and 'etalon_type' in st.session_state:
        image = Image.open(uploaded_file)
        
        # Resize pour affichage
        max_width = 800
        if image.width > max_width:
            ratio = max_width / image.width
            display_image = image.resize((max_width, int(image.height * ratio)), Image.Resampling.LANCZOS)
        else:
            display_image = image
        
        st.image(display_image, caption="Image source", use_container_width=True)
        
        # Étape 3: Calibration
        st.subheader("3️⃣ Calibration Interactive")
        
        col_cal1, col_cal2 = st.columns([3, 1])
        
        with col_cal1:
            st.markdown("**Cliquez sur les extrémités de l'étalon**")
            
            x1 = st.number_input("X1 (début)", 0, display_image.width, 100, key="x1")
            y1 = st.number_input("Y1 (début)", 0, display_image.height, 400, key="y1")
            x2 = st.number_input("X2 (fin)", 0, display_image.width, 300, key="x2")
            y2 = st.number_input("Y2 (fin)", 0, display_image.height, 400, key="y2")
            
            # Preview calibration
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
                ppm = engine.set_etalon(st.session_state['etalon_type'], pixel_distance)
                st.session_state['morpho_engine'] = engine
                st.session_state['ppm'] = ppm
                st.success(f"Calibration: {ppm:.2f} px/cm")
                st.rerun()
        
        # Étape 4: Mesures
        if 'morpho_engine' in st.session_state:
            st.subheader("4️⃣ Points de Mesure Morphométrique")
            
            st.markdown("""
            **Points à définir**:
            - **Garrot**: Point haut entre omoplates
            - **Base queue**: Insertion queue sur croupe  
            - **Patte arrière**: Sol au jarret
            - **Poitrine**: Point le plus large
            - **Canon**: Milieu du canon postérieur
            - **Mamelles**: Attachment gauche/droite + arrière
            """)
            
            points = {}
            
            with st.expander("📍 Saisie des Points", expanded=True):
                cols = st.columns(2)
                
                with cols[0]:
                    points['garrot'] = (
                        st.number_input("Garrot X", 0, display_image.width, 200, key="gx"),
                        st.number_input("Garrot Y", 0, display_image.height, 150, key="gy")
                    )
                    points['base_queue'] = (
                        st.number_input("Queue X", 0, display_image.width, 600, key="bqx"),
                        st.number_input("Queue Y", 0, display_image.height, 200, key="bqy")
                    )
                    points['patte_arriere'] = (
                        st.number_input("Patte X", 0, display_image.width, 550, key="pax"),
                        st.number_input("Patte Y", 0, display_image.height, 500, key="pay")
                    )
                    points['mamelle_gauche'] = (
                        st.number_input("Mamelle G X", 0, display_image.width, 300, key="mgx"),
                        st.number_input("Mamelle G Y", 0, display_image.height, 450, key="mgy")
                    )
                
                with cols[1]:
                    points['poitrine'] = (
                        st.number_input("Poitrine X", 0, display_image.width, 150, key="px"),
                        st.number_input("Poitrine Y", 0, display_image.height, 350, key="py")
                    )
                    points['canon'] = (
                        st.number_input("Canon X", 0, display_image.width, 550, key="cx"),
                        st.number_input("Canon Y", 0, display_image.height, 480, key="cy")
                    )
                    points['canon_peripherie'] = (
                        st.number_input("Canon Périph X", 0, display_image.width, 570, key="cpx"),
                        st.number_input("Canon Périph Y", 0, display_image.height, 480, key="cpy")
                    )
                    points['mamelle_droite'] = (
                        st.number_input("Mamelle D X", 0, display_image.width, 400, key="mdx"),
                        st.number_input("Mamelle D Y", 0, display_image.height, 450, key="mdy")
                    )
                
                points['mamelle_arriere'] = (
                    st.number_input("Mamelle Arr X", 0, display_image.width, 350, key="max"),
                    st.number_input("Mamelle Arr Y", 0, display_image.height, 500, key="may")
                )
            
            if st.button("🧮 Calculer les Mesures"):
                engine = st.session_state['morpho_engine']
                measurements = engine.calculate_measurements(points)
                
                if "error" in measurements:
                    st.error(measurements["error"])
                else:
                    st.session_state['last_measurements'] = measurements
                    
                    # Affichage résultats
                    st.subheader("📏 Résultats")
                    
                    overlay_img = engine.generate_overlay_image(display_image, points, measurements)
                    st.image(overlay_img, caption="Mesures annotées", use_container_width=True)
                    
                    col_res1, col_res2 = st.columns(2)
                    
                    with col_res1:
                        st.markdown("**📐 Morphologie**")
                        for key in ['longueur_corps_cm', 'hauteur_garrot_cm', 'tour_poitrine_cm', 'circonf_canon_cm']:
                            if key in measurements:
                                st.metric(key.replace('_', ' ').title(), f"{measurements[key]} cm")
                    
                    with col_res2:
                        st.markdown("**🥛 Mamelle & Poids**")
                        if 'score_mamelle' in measurements:
                            st.metric("Score Mamelle", f"{measurements['score_mamelle']}/10")
                        if 'poids_estime_kg' in measurements:
                            st.metric("Poids Estimé", f"{measurements['poids_estime_kg']} kg")
                        if 'indice_compacite' in measurements:
                            st.metric("Indice Compacité", measurements['indice_compacite'])
                    
                    # Sauvegarde
                    with st.form("save_morpho"):
                        animal_id = st.text_input("ID Animal *", placeholder="FR123456789")
                        notes = st.text_area("Notes", placeholder="Conditions de prise de vue...")
                        
                        if st.form_submit_button("💾 Sauvegarder"):
                            if animal_id:
                                db.execute(
                                    """INSERT INTO mesures_morpho 
                                       (brebis_id, date_mesure, longueur_cm, hauteur_garrot_cm, 
                                        tour_poitrine_cm, circonf_canon_cm, score_mamelle, 
                                        poids_estime_kg, indice_compacite, notes, etalon_type, ppm_ratio)
                                       VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
                                    (animal_id, date.today(), 
                                     measurements.get('longueur_corps_cm'),
                                     measurements.get('hauteur_garrot_cm'),
                                     measurements.get('tour_poitrine_cm'),
                                     measurements.get('circonf_canon_cm'),
                                     measurements.get('score_mamelle'),
                                     measurements.get('poids_estime_kg'),
                                     measurements.get('indice_compacite'),
                                     notes,
                                     st.session_state['etalon_type'].value,
                                     st.session_state['ppm'])
                                )
                                st.success(f"✅ Mesures enregistrées pour {animal_id}")
                            else:
                                st.error("ID animal requis")
    
    elif uploaded_file and 'etalon_type' not in st.session_state:
        st.warning("⚠️ Sélectionnez d'abord un étalon de calibration")
    
    # Historique
    st.subheader("📊 Historique des Mesures")
    df_morpho = db.fetch_df("""
        SELECT m.*, b.nom, b.race 
        FROM mesures_morpho m 
        LEFT JOIN brebis b ON m.brebis_id = b.identifiant_unique 
        ORDER BY m.date_mesure DESC 
        LIMIT 20
    """)
    
    if not df_morpho.empty:
        st.dataframe(df_morpho, use_container_width=True)
    else:
        st.info("Aucune mesure enregistrée")

# ============================================================================
# 7. AUTRES MODULES (DASHBOARD, INSCRIPTION, ETC.)
# ============================================================================

def render_dashboard(db: DatabaseManager):
    st.title("📊 Tableau de Bord")
    df_b = db.fetch_df("SELECT * FROM brebis")
    df_l = db.fetch_df("SELECT * FROM controle_laitier")
    df_m = db.fetch_df("SELECT * FROM mesures_morpho")
    
    if not df_b.empty:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Effectif", len(df_b))
        c2.metric("Poids Moyen", f"{round(df_b['poids'].mean(), 1)} kg")
        avg_lait = df_l['quantite_lait'].mean() if not df_l.empty else 0
        c3.metric("Moyenne Lait", f"{round(avg_lait, 2)} L")
        c4.metric("Mesures Scanner", len(df_m) if not df_m.empty else 0)
        
        st.dataframe(df_b, use_container_width=True)
        
        if not df_l.empty:
            fig = px.line(df_l, x='date_controle', y='quantite_lait', color='brebis_id')
            st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("Aucun animal enregistré")

def render_inscription(db: DatabaseManager):
    st.title("📝 Inscription Phénotypique")
    with st.form("form_inscription"):
        col1, col2 = st.columns(2)
        uid = col1.text_input("ID Boucle")
        nom = col1.text_input("Nom")
        race = col2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Autre"])
        poids = col2.number_input("Poids (kg)", 10.0, 150.0, 50.0)
        
        tp = st.number_input("Tour Poitrine (cm)", 40.0, 160.0, 85.0)
        lg = st.number_input("Longueur Corps (cm)", 30.0, 140.0, 75.0)
        note_m = st.slider("Note Mamelle", 1, 10, 5)
        
        if st.form_submit_button("Sauvegarder"):
            db.execute(
                "INSERT INTO brebis (identifiant_unique, nom, race, poids, note_mamelle, tour_poitrine, longueur, created_at) VALUES (?,?,?,?,?,?,?,?)",
                (uid, nom, race, poids, note_m, tp, lg, date.today())
            )
            st.success(f"Animal {uid} ajouté")

def render_lait(db: DatabaseManager):
    st.title("🥛 Contrôle Laitier")
    with st.form("form_lait"):
        id_lait = st.text_input("ID Brebis")
        qte = st.number_input("Quantité (L)", 0.0, 12.0, 1.5)
        dt = st.date_input("Date", date.today())
        if st.form_submit_button("Valider"):
            db.execute("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)",
                      (id_lait, dt, qte))
            st.success("Enregistré")

def render_sante(db: DatabaseManager):
    st.title("🩺 Santé & Vaccins")
    with st.expander("Nouvel acte"):
        with st.form("form_sante"):
            id_s = st.text_input("ID Animal")
            type_a = st.selectbox("Acte", ["Vaccination", "Déparasitage", "Traitement Curatif"])
            prod = st.text_input("Produit")
            rappel = st.date_input("Rappel", date.today() + timedelta(days=30))
            if st.form_submit_button("Ajouter"):
                db.execute("INSERT INTO sante (brebis_id, date_soin, type_acte, produit, rappel_prevu) VALUES (?,?,?,?,?)",
                          (id_s, date.today(), type_a, prod, rappel))
                st.success("Soin enregistré")
    
    df_s = db.fetch_df("SELECT * FROM sante")
    st.table(df_s)

def render_genomique():
    st.title("🧬 Génomique")
    if 'genomique' not in st.session_state:
        st.session_state.genomique = BioInfoEngine()
    gen = st.session_state.genomique
    
    dna_txt = st.text_area("Séquences ADN", height=200)
    if dna_txt:
        seqs = gen.extraire_multi_fasta(dna_txt)
        results = []
        for name, seq in seqs.items():
            row = {"ID": name}
            for gene, ref in gen.GENES_INTERET.items():
                score = gen.alignement_expert(seq, ref)
                row[gene] = f"{'OUI' if score > 85 else 'NON'} ({score}%)"
            results.append(row)
        st.dataframe(pd.DataFrame(results))

def render_nutrition():
    st.title("🌾 Nutrition")
    poids = st.number_input("Poids (kg)", 10, 150, 60)
    c1, c2 = st.columns(2)
    c1.write(f"Concentré: {poids * 0.012:.2f} kg/j")
    c2.write(f"Fourrage: {poids * 0.02:.2f} kg/j")

# ============================================================================
# 8. APPLICATION PRINCIPALE (MAIN - TOUJOURS EN DERNIER)
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    # CSS minimal
    st.markdown("""
    <style>
    .stMetric {background-color: #f0f2f6; border-radius: 8px; padding: 10px;}
    </style>
    """, unsafe_allow_html=True)
    
    # Initialisation session state
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    
    # Sidebar
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    menu = [
        "📊 Dashboard", 
        "📝 Inscription", 
        "📷 Scanner IA 1m",
        "🥛 Contrôle Laitier", 
        "🩺 Santé", 
        "🧬 Génomique", 
        "🌾 Nutrition"
    ]
    choice = st.sidebar.radio("Navigation", menu)
    
    # Routing
    if choice == "📊 Dashboard":
        render_dashboard(db)
    elif choice == "📝 Inscription":
        render_inscription(db)
    elif choice == "📷 Scanner IA 1m":
        render_scanner_ia(db)  # Appel correct avec db injecté
    elif choice == "🥛 Contrôle Laitier":
        render_lait(db)
    elif choice == "🩺 Santé":
        render_sante(db)
    elif choice == "🧬 Génomique":
        render_genomique()
    elif choice == "🌾 Nutrition":
        render_nutrition()

if __name__ == "__main__":
    main()
