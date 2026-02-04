"""
EXPERT OVIN DZ PRO - Version Ultime Complète
Avec:
- Export vers R (rpy2)
- Templates PDF protocoles terrain
- Synchronisation cloud (Firebase/AWS)
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
import tempfile
import os
import hashlib
import uuid
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, asdict
from enum import Enum
from pathlib import Path

# Tentative import rpy2 pour export R
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False
    st.warning("rpy2 non installé - Export R désactivé (utiliser CSV/Excel fallback)")

# Imports pour PDF
try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image as RLImage
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.pdfgen import canvas
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    REPORTLAB_AVAILABLE = True
except ImportError:
    REPORTLAB_AVAILABLE = False

# Imports Firebase/AWS
try:
    import firebase_admin
    from firebase_admin import credentials, firestore, storage
    FIREBASE_AVAILABLE = True
except ImportError:
    FIREBASE_AVAILABLE = False

try:
    import boto3
    from botocore.exceptions import ClientError
    AWS_AVAILABLE = False  # Activé manuellement si credentials configurés
except ImportError:
    AWS_AVAILABLE = False

# OpenCV et traitement d'image
try:
    import cv2
    import numpy as np_cv
    from PIL import Image, ImageDraw, ImageFont
    OPENCV_AVAILABLE = True
except ImportError:
    OPENCV_AVAILABLE = False

# Scikit-image
try:
    from skimage import segmentation, filters, measure, morphology
    from skimage.color import rgb2gray, rgb2hsv
    from skimage.feature import canny
    from scipy import ndimage, stats as sci_stats
    SKIMAGE_AVAILABLE = True
except ImportError:
    SKIMAGE_AVAILABLE = False

# Statsmodels
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.diagnostic import het_breuschpagan, normal_ad
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False

# Scikit-learn
try:
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.ensemble import RandomForestRegressor
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# ============================================================================
# SECTION 2: CONFIGURATION ET CONSTANTES
# ============================================================================
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constantes morphométriques
REFERENCE_MORPHO_OVINE = {
    'hauteur_garrot_cm': {'min': 60, 'max': 75, 'optimal': 68, 'unit': 'cm'},
    'longueur_corps_cm': {'min': 70, 'max': 85, 'optimal': 78, 'unit': 'cm'},
    'largeur_bassin_cm': {'min': 18, 'max': 24, 'optimal': 21, 'unit': 'cm'},
    'tour_poitrine_cm': {'min': 85, 'max': 100, 'optimal': 92, 'unit': 'cm'},
    'longueur_mamelle_cm': {'min': 15, 'max': 22, 'optimal': 18, 'unit': 'cm'},
    'largeur_mamelle_cm': {'min': 12, 'max': 18, 'optimal': 15, 'unit': 'cm'},
    'profondeur_mamelle_cm': {'min': 12, 'max': 20, 'optimal': 16, 'unit': 'cm'},
    'ecart_tetines_cm': {'min': 8, 'max': 14, 'optimal': 11, 'unit': 'cm'},
}

# ============================================================================
# MODULE: EXPORT VERS R (RPY2)
# ============================================================================

class RExportManager:
    """
    Gestion export données et analyses vers R véritable
    Utilise rpy2 pour interopérabilité Python-R
    """
    
    def __init__(self):
        self.r_available = RPY2_AVAILABLE
        self.r_libs_loaded = False
        
    def init_r_environment(self):
        """Initialise l'environnement R avec packages nécessaires"""
        if not self.r_available:
            return False
        
        try:
            # Activation conversion pandas-R
            pandas2ri.activate()
            
            # Import packages R de base
            self.r_base = importr('base')
            self.r_stats = importr('stats')
            self.r_utils = importr('utils')
            
            # Tentative chargement packages additionnels
            try:
                self.r_ggplot2 = importr('ggplot2')
                self.r_dplyr = importr('dplyr')
                self.r_corrplot = importr('corrplot')
                self.r_facto = importr('FactoMineR')
                self.r_libs_loaded = True
            except Exception as e:
                st.warning(f"Packages R optionnels non chargés: {e}")
                self.r_libs_loaded = True  # Base suffisant pour export
                
            return True
            
        except Exception as e:
            st.error(f"Erreur initialisation R: {e}")
            return False
    
    def export_dataframe_to_r(self, df: pd.DataFrame, r_varname: str = "df_ovins") -> bool:
        """Exporte DataFrame pandas vers variable R"""
        if not self.r_available:
            return False
        
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_df = ro.conversion.py2rpy(df)
                ro.globalenv[r_varname] = r_df
            return True
        except Exception as e:
            st.error(f"Erreur export vers R: {e}")
            return False
    
    def generate_r_script(self, df: pd.DataFrame, analysis_type: str = "complete") -> str:
        """
        Génère script R complet prêt à exécuter
        Alternative si rpy2 non disponible
        """
        
        script = f"""# Script R généré automatiquement par Expert Ovin DZ Pro
# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# Données: {len(df)} observations, {len(df.columns)} variables

# ============================================================================
# 1. CHARGEMENT ET PRÉPARATION
# ============================================================================

# Lecture données (remplacer chemin si nécessaire)
# df <- read.csv("votre_fichier.csv", stringsAsFactors = FALSE)
# Ou si déjà en mémoire:
# df <- df_ovins

# Structure données
str(df)
summary(df)

# Vérification valeurs manquantes
library(visdat)
vis_miss(df)

# ============================================================================
# 2. STATISTIQUES DESCRIPTIVES
# ============================================================================

library(summarytools)
dfSummary(df)

# Stats spécifiques morphométrie
library(pastecs)
stat.desc(df[, sapply(df, is.numeric)], norm = TRUE)

# ============================================================================
# 3. CORRÉLATIONS
# ============================================================================

# Matrice corrélation
M <- cor(df[, sapply(df, is.numeric)], use = "complete.obs")

# Visualisation
library(corrplot)
corrplot(M, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         diag = FALSE)

# Test significativité
library(Hmisc)
res <- rcorr(as.matrix(df[, sapply(df, is.numeric)]))
print(res$P)  # p-values

# ============================================================================
# 4. RÉGRESSION - Prédiction Production Laitière
# ============================================================================

# Modèle linéaire multiple
model_lm <- lm(production_lait_jour ~ longueur_mamelle + profondeur_mamelle + 
               score_symetrie + score_attache_mamelle, data = df)

summary(model_lm)

# Diagnostics
par(mfrow = c(2, 2))
plot(model_lm)

# Test multicolinéarité
library(car)
vif(model_lm)

# Régression robuste si hétéroscédasticité
library(MASS)
model_rlm <- rlm(production_lait_jour ~ ., data = df)

# ============================================================================
# 5. ANALYSE MULTIVARIÉE - PCA
# ============================================================================

library(FactoMineR)
library(factoextra)

# Sélection variables morphométriques
vars_morpho <- c("hauteur_garrot", "longueur_corps", "largeur_bassin",
                 "longueur_mamelle", "largeur_mamelle", "profondeur_mamelle")

res_pca <- PCA(df[, vars_morpho], scale.unit = TRUE, graph = FALSE)

# Visualisations
fviz_eig(res_pca, addlabels = TRUE)  # Scree plot
fviz_pca_ind(res_pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_pca_var(res_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# ============================================================================
# 6. CLUSTERING - Classification morphologique
# ============================================================================

# Préparation données
df_scaled <- scale(df[, vars_morpho])

# Détermination nombre optimal clusters (méthode elbow)
library(factoextra)
fviz_nbclust(df_scaled, kmeans, method = "wss")

# K-means avec k=3 (exemple)
set.seed(123)
km_res <- kmeans(df_scaled, centers = 3, nstart = 25)

# Visualisation
fviz_cluster(km_res, data = df_scaled, palette = "jco",
             ggtheme = theme_minimal())

# Profil clusters
df$cluster <- as.factor(km_res$cluster)
aggregate(df[, vars_morpho], by = list(cluster = df$cluster), FUN = mean)

# ============================================================================
# 7. MODÈLE MIXTE - Si données répétées
# ============================================================================

library(lme4)
library(lmerTest)

# Exemple: effet aléatoire animal sur mesures répétées
# model_mixte <- lmer(production_lait_jour ~ (1|animal_id) + age_mois, data = df)

# ============================================================================
# 8. EXPORT RÉSULTATS
# ============================================================================

# Sauvegarde environnement
save.image(file = "session_expert_ovin.RData")

# Export résultats principaux
write.csv(summary(model_lm)$coefficients, "coefficients_regression.csv")
write.csv(M, "matrice_correlation.csv")

# Génération rapport automatique
library(rmarkdown)
# render("rapport_ovin.Rmd", output_format = "pdf_document")

print("Analyse terminée avec succès!")
"""
        return script
    
    def execute_r_analysis(self, df: pd.DataFrame, analysis_type: str = "descriptive") -> Dict:
        """
        Exécute directement analyses R via rpy2
        Retourne résultats sous forme de dictionnaire
        """
        if not self.init_r_environment():
            return {"error": "R non disponible"}
        
        results = {}
        
        try:
            # Export données vers R
            self.export_dataframe_to_r(df, "df_ovins")
            
            if analysis_type == "descriptive":
                # Stats descriptives
                r_code = """
                summary_stats <- summary(df_ovins)
                summary_stats
                """
                results['summary'] = ro.r(r_code)
                
            elif analysis_type == "correlation":
                # Matrice corrélation
                r_code = """
                cor_matrix <- cor(df_ovins[, sapply(df_ovins, is.numeric)], use = "complete.obs")
                cor_matrix
                """
                results['correlation'] = ro.r(r_code)
                
            elif analysis_type == "regression":
                # Régression simple exemple
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) >= 2:
                    y = numeric_cols[0]
                    x = numeric_cols[1]
                    r_code = f"""
                    model <- lm({y} ~ {x}, data = df_ovins)
                    summary(model)
                    """
                    results['regression'] = ro.r(r_code)
            
            return results
            
        except Exception as e:
            return {"error": str(e)}
    
    def render_r_export_interface(self, df: pd.DataFrame):
        """Interface Streamlit pour export R"""
        st.markdown("## 🧮 Export vers R (Statistiques Avancées)")
        
        if not self.r_available:
            st.warning("""
            **rpy2 non installé** - Mode fallback CSV/R-script activé
            
            Installation: `pip install rpy2` (requiert R installé système)
            """)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📤 Export Données")
            
            if self.r_available:
                if st.button("Envoyer vers R (mémoire)"):
                    if self.export_dataframe_to_r(df, "df_ovins"):
                        st.success("✅ DataFrame transféré vers R (variable: df_ovins)")
                        
                        # Affichage preview R
                        preview = ro.r("head(df_ovins, 3)")
                        st.text("Aperçu dans R:")
                        st.text(str(preview))
            else:
                st.info("Mode CSV activé")
            
            # Export CSV universel
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            st.download_button(
                "📥 Télécharger CSV (pour R)",
                csv_buffer.getvalue(),
                f"expert_ovin_data_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv"
            )
        
        with col2:
            st.subheader("📜 Script R Complet")
            
            script_type = st.selectbox("Type d'analyse", 
                                      ["Complète", "Descriptive", "Régression", "PCA", "Clustering"])
            
            script = self.generate_r_script(df, script_type.lower())
            
            st.download_button(
                "📥 Télécharger Script .R",
                script,
                f"analyse_ovin_{script_type.lower()}.R",
                mime="text/plain"
            )
            
            with st.expander("Voir le script"):
                st.code(script, language='r')
        
        # Exécution directe si R disponible
        if self.r_available and st.checkbox("Exécuter analyse R directement"):
            analysis = st.selectbox("Analyse à exécuter", 
                                   ["Descriptive", "Corrélation", "Régression"])
            
            if st.button("🚀 Lancer dans R"):
                with st.spinner("Exécution R..."):
                    results = self.execute_r_analysis(df, analysis.lower())
                    
                    if 'error' in results:
                        st.error(results['error'])
                    else:
                        st.success("Analyse R terminée")
                        for key, val in results.items():
                            st.subheader(key.capitalize())
                            st.text(str(val))

# ============================================================================
# MODULE: TEMPLATES PDF TERRAIN
# ============================================================================

class PDFTemplateManager:
    """
    Génération protocoles PDF imprimables pour terrain
    Fiches de mesure, guides opérateur, rapports
    """
    
    def __init__(self):
        self.available = REPORTLAB_AVAILABLE
        
    def generate_fiche_mesure_traditionnelle(self) -> bytes:
        """
        Génère fiche A4 traditionnelle avec cases à cocher/cocher
        Pour utilisation sur le terrain sans tablette
        """
        if not self.available:
            return b""
        
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=A4,
                               rightMargin=2*cm, leftMargin=2*cm,
                               topMargin=2*cm, bottomMargin=2*cm)
        
        elements = []
        styles = getSampleStyleSheet()
        
        # Style titre
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=18,
            textColor=colors.HexColor('#2E5090'),
            spaceAfter=30,
            alignment=TA_CENTER
        )
        
        # En-tête
        elements.append(Paragraph("EXPERT OVIN DZ PRO", title_style))
        elements.append(Paragraph("<b>FICHE DE MESURE MORPHOMÉTRIQUE</b>", styles['Heading2']))
        elements.append(Spacer(1, 0.5*cm))
        
        # Informations générales
        data_info = [
            ["Date:", "____/____/______", "Heure:", "____:____"],
            ["Opérateur:", "_________________________", "ID Animal:", "_________________________"],
            ["Race:", "_________________________", "Âge (mois):", "_______"],
            ["N° Lactation:", "_______", "Élevage:", "_________________________"]
        ]
        
        table_info = Table(data_info, colWidths=[4*cm, 6*cm, 4*cm, 4*cm])
        table_info.setStyle(TableStyle([
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ]))
        elements.append(table_info)
        elements.append(Spacer(1, 1*cm))
        
        # Mensurations corporelles
        elements.append(Paragraph("<b>📐 MENSURATIONS CORPORELLES (au cm près)</b>", styles['Heading3']))
        
        data_corp = [
            ["Mesure", "Valeur", "Fourchette réf.", "Observation"],
            ["Hauteur au garrot", "_______ cm", "60-75 cm", "□"],
            ["Longueur du corps", "_______ cm", "70-85 cm", "□"],
            ["Largeur du bassin", "_______ cm", "18-24 cm", "□"],
            ["Tour de poitrine", "_______ cm", "85-100 cm", "□"],
            ["Profondeur poitrine", "_______ cm", "30-38 cm", "□"],
            ["Largeur épaules", "_______ cm", "-", "□"],
        ]
        
        table_corp = Table(data_corp, colWidths=[6*cm, 3*cm, 4*cm, 3*cm])
        table_corp.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ]))
        elements.append(table_corp)
        elements.append(Spacer(1, 0.8*cm))
        
        # Mensurations mamelle
        elements.append(Paragraph("<b>🥛 MENSURATIONS MAMELLE (au cm près)</b>", styles['Heading3']))
        
        data_mam = [
            ["Mesure", "Valeur", "Fourchette réf.", "Observation"],
            ["Longueur mamelle", "_______ cm", "15-22 cm", "□"],
            ["Largeur mamelle", "_______ cm", "12-18 cm", "□"],
            ["Profondeur mamelle", "_______ cm", "12-20 cm", "□"],
            ["Circonf. tétine G", "_______ cm", "2.5-4.0 cm", "□"],
            ["Circonf. tétine D", "_______ cm", "2.5-4.0 cm", "□"],
            ["Écart tétines", "_______ cm", "8-14 cm", "□"],
            ["Hauteur attache", "_______ cm", "-", "□"],
            ["Angle attache", "_______ °", "30-90°", "□"],
            ["Angle tétines", "_______ °", "0-45°", "□"],
        ]
        
        table_mam = Table(data_mam, colWidths=[6*cm, 3*cm, 4*cm, 3*cm])
        table_mam.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#E06666')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#FCE4D6')),
        ]))
        elements.append(table_mam)
        elements.append(Spacer(1, 0.8*cm))
        
        # Scores subjectifs
        elements.append(Paragraph("<b>⭐ SCORES SUBJECTIFS (échelle 1-9)</b>", styles['Heading3']))
        
        data_scores = [
            ["Critère", "Score (1-9)", "Description"],
            ["Attache mamelle", "_____", "1=Lâche → 9=Serrée"],
            ["Profondeur mamelle", "_____", "1=Superficielle → 9=Profonde"],
            ["Symétrie", "_____", "1=Asymétrique → 9=Parfaite"],
            ["État corporel", "_____", "1=Maigre → 5=Optimal → 9=Obèse"],
        ]
        
        table_scores = Table(data_scores, colWidths=[6*cm, 3*cm, 8*cm])
        table_scores.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#70AD47')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ]))
        elements.append(table_scores)
        elements.append(Spacer(1, 1*cm))
        
        # Section production
        elements.append(Paragraph("<b>🥛 PRODUCTION LAITIÈRE (si disponible)</b>", styles['Heading3']))
        
        data_prod = [
            ["Production jour", "_______ L", "Date contrôle", "____/____/______"],
            ["% Matière grasse", "_______ %", "% Protéines", "_______ %"],
        ]
        
        table_prod = Table(data_prod, colWidths=[6*cm, 3*cm, 6*cm, 4*cm])
        table_prod.setStyle(TableStyle([
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
        ]))
        elements.append(table_prod)
        elements.append(Spacer(1, 1.5*cm))
        
        # Pied de page
        elements.append(Paragraph("<i>Protocole ISO 7481 - Expert Ovin DZ Pro</i>", styles['Italic']))
        elements.append(Paragraph(f"<i>Généré le {datetime.now().strftime('%d/%m/%Y')}</i>", styles['Italic']))
        
        # Génération PDF
        doc.build(elements)
        pdf = buffer.getvalue()
        buffer.close()
        return pdf
    
    def generate_guide_operateur(self) -> bytes:
        """
        Guide méthodologique PDF pour opérateurs terrain
        """
        if not self.available:
            return b""
        
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=A4)
        elements = []
        styles = getSampleStyleSheet()
        
        # Contenu pédagogique
        title = Paragraph("<b>GUIDE DE L'OPÉRATEUR - PROTOCOLE MESURE</b>", 
                         ParagraphStyle('Title', fontSize=20, textColor=colors.HexColor('#1F4E78'), 
                                       alignment=TA_CENTER, spaceAfter=20))
        elements.append(title)
        
        # Sections avec illustrations textuelles
        sections = [
            ("1. POSITIONNEMENT DE L'ANIMAL", 
             "L'animal doit être en station debout naturelle, quatre membres bien répartis, "
             "sur sol plat et dur. Attendre la fin des mouvements de mastication."),
            ("2. HAUTEUR AU GARROT", 
             "Mesurer verticalement du point le plus haut du dos (dernière vertèbre thoracique) "
             "jusqu'au sol. Utiliser un bâton gradué ou mesure laser."),
            ("3. MAMELLE - ATTENTION PARTICULIÈRE",
             "Ne pas comprimer le tissu mammaire. Mesurer avant la traite du matin. "
             "Noter toute asymétrie visible ou anomalie de forme."),
        ]
        
        for titre, contenu in sections:
            elements.append(Paragraph(f"<b>{titre}</b>", styles['Heading3']))
            elements.append(Paragraph(contenu, styles['BodyText']))
            elements.append(Spacer(1, 0.5*cm))
        
        # Tableau erreurs fréquentes
        data_erreurs = [
            ["Erreur fréquente", "Conséquence", "Correction"],
            ["Animal mal positionné", "Biais +2-3 cm", "Attendre station naturelle"],
            ["Ruban trop serré", "Sous-estimation volume", "Contact sans compression"],
            ["Mesure post-traite", "Volume sous-estimé", "Mesurer avant traite"],
        ]
        
        table_err = Table(data_erreurs, colWidths=[6*cm, 6*cm, 6*cm])
        table_err.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#C00000')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ]))
        elements.append(table_err)
        
        doc.build(elements)
        return buffer.getvalue()
    
    def generate_rapport_pdf(self, mesure: 'MesuresManuelles', indices: Dict) -> bytes:
        """
        Génère rapport PDF individuel pour un animal
        """
        if not self.available:
            return b""
        
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=A4)
        elements = []
        styles = getSampleStyleSheet()
        
        # Rapport stylisé avec résultats
        elements.append(Paragraph(f"<b>RAPPORT MORPHOMÉTRIQUE</b>", 
                                 ParagraphStyle('Title', fontSize=18, alignment=TA_CENTER)))
        elements.append(Paragraph(f"Animal: <b>{mesure.animal_id}</b> | Race: {mesure.race}", 
                                 styles['Heading2']))
        
        # Tableau récapitulatif
        data_recap = [
            ["Paramètre", "Mesuré", "Référence", "Statut"],
            ["Hauteur garrot", f"{mesure.hauteur_garrot} cm", "68 cm", 
             "✓" if 60 <= mesure.hauteur_garrot <= 75 else "⚠"],
            ["Longueur corps", f"{mesure.longueur_corps} cm", "78 cm",
             "✓" if 70 <= mesure.longueur_corps <= 85 else "⚠"],
            ["Volume mamelle", f"{indices.get('volume_mamelle_cm3', 0)} cm³", "-", "ℹ"],
            ["Score global", f"{indices.get('score_morphologique_global', 0)}/9", ">7", 
             "✓" if indices.get('score_morphologique_global', 0) >= 7 else "⚠"],
        ]
        
        table_recap = Table(data_recap, colWidths=[5*cm, 3*cm, 3*cm, 2*cm])
        table_recap.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4472C4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
        ]))
        elements.append(table_recap)
        
        # Interprétation
        elements.append(Spacer(1, 1*cm))
        elements.append(Paragraph("<b>Interprétation:</b>", styles['Heading3']))
        
        score = indices.get('score_morphologique_global', 0)
        if score >= 8:
            interpretation = "Excellente morphologie mamelle - Recommandé pour reproduction"
            color = colors.green
        elif score >= 6:
            interpretation = "Bonne morphologie - Surveillance recommandée"
            color = colors.orange
        else:
            interpretation = "Morphologie à améliorer - Réforme à considérer"
            color = colors.red
        
        elements.append(Paragraph(interpretation, 
                                 ParagraphStyle('Interp', textColor=color, fontSize=12)))
        
        doc.build(elements)
        return buffer.getvalue()
    
    def render_pdf_interface(self):
        """Interface Streamlit pour génération PDF"""
        st.markdown("## 📄 Générateur de Documents PDF")
        
        if not self.available:
            st.error("ReportLab requis: pip install reportlab")
            return
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("Fiche de Mesure")
            st.caption("Formulaire vierge pour terrain")
            if st.button("Générer fiche A4"):
                pdf = self.generate_fiche_mesure_traditionnelle()
                st.download_button("📥 Télécharger PDF", pdf, 
                                 f"fiche_mesure_{datetime.now().strftime('%Y%m%d')}.pdf",
                                 mime="application/pdf")
        
        with col2:
            st.subheader("Guide Opérateur")
            st.caption("Protocole méthodologique")
            if st.button("Générer guide"):
                pdf = self.generate_guide_operateur()
                st.download_button("📥 Télécharger PDF", pdf,
                                 "guide_operateur_ovin.pdf",
                                 mime="application/pdf")
        
        with col3:
            st.subheader("Rapport Individuel")
            st.caption("Nécessite données saisies")
            
            # Vérifier données disponibles
            mesures = st.session_state.get('mesures_manuelles', [])
            if mesures:
                options = [f"{i+1}. {m['mesure'].animal_id}" for i, m in enumerate(mesures)]
                selection = st.selectbox("Sélection animal", options)
                
                if st.button("Générer rapport"):
                    idx = int(selection.split('.')[0]) - 1
                    mesure = mesures[idx]['mesure']
                    indices = mesures[idx]['indices']
                    pdf = self.generate_rapport_pdf(mesure, indices)
                    st.download_button("📥 Télécharger rapport", pdf,
                                     f"rapport_{mesure.animal_id}.pdf",
                                     mime="application/pdf")
            else:
                st.info("Aucune mesure disponible")

# ============================================================================
# MODULE: SYNCHRONISATION CLOUD
# ============================================================================

class CloudSyncManager:
    """
    Gestion synchronisation données cloud
    Support Firebase (Firestore + Storage) et AWS S3
    """
    
    def __init__(self):
        self.firebase_initialized = False
        self.aws_client = None
        self.local_db_path = "expert_ovin_local.db"
        self._init_local_db()
        
    def _init_local_db(self):
        """Initialise SQLite local pour cache offline"""
        conn = sqlite3.connect(self.local_db_path)
        cursor = conn.cursor()
        
        # Table mesures
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS mesures_morpho (
                id TEXT PRIMARY KEY,
                animal_id TEXT,
                date_mesure TEXT,
                data_json TEXT,
                synced INTEGER DEFAULT 0,
                cloud_id TEXT
            )
        ''')
        
        # Table images
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS images_analyse (
                id TEXT PRIMARY KEY,
                animal_id TEXT,
                image_blob BLOB,
                date_capture TEXT,
                synced INTEGER DEFAULT 0
            )
        ''')
        
        conn.commit()
        conn.close()
    
    def init_firebase(self, credentials_dict: Dict = None):
        """Initialise connexion Firebase"""
        if not FIREBASE_AVAILABLE:
            st.error("Firebase non installé: pip install firebase-admin")
            return False
        
        try:
            if credentials_dict:
                # Mode credentials dict (streamlit secrets)
                cred = credentials.Certificate(credentials_dict)
            else:
                # Mode fichier
                cred = credentials.Certificate("firebase-credentials.json")
            
            firebase_admin.initialize_app(cred, {
                'storageBucket': 'your-bucket.appspot.com'
            })
            self.firebase_initialized = True
            self.db = firestore.client()
            self.bucket = storage.bucket()
            return True
            
        except Exception as e:
            st.error(f"Erreur Firebase: {e}")
            return False
    
    def init_aws(self, aws_access_key: str, aws_secret_key: str, region: str = "eu-west-3"):
        """Initialise connexion AWS S3"""
        if not AWS_AVAILABLE:
            st.error("Boto3 non installé: pip install boto3")
            return False
        
        try:
            self.aws_client = boto3.client(
                's3',
                aws_access_key_id=aws_access_key,
                aws_secret_access_key=aws_secret_key,
                region_name=region
            )
            # Test connexion
            self.aws_client.list_buckets()
            return True
        except Exception as e:
            st.error(f"Erreur AWS: {e}")
            return False
    
    def save_measurement_local(self, mesure: 'MesuresManuelles', indices: Dict) -> str:
        """Sauvegarde locale SQLite (toujours disponible)"""
        conn = sqlite3.connect(self.local_db_path)
        cursor = conn.cursor()
        
        record_id = str(uuid.uuid4())
        data = {
            'mesure': mesure.to_dict(),
            'indices': indices,
            'timestamp': datetime.now().isoformat()
        }
        
        cursor.execute('''
            INSERT INTO mesures_morpho (id, animal_id, date_mesure, data_json, synced)
            VALUES (?, ?, ?, ?, 0)
        ''', (record_id, mesure.animal_id, mesure.date_mesure.isoformat(), 
              json.dumps(data)))
        
        conn.commit()
        conn.close()
        return record_id
    
    def sync_to_firebase(self, record_id: str = None) -> Dict:
        """Synchronise données locales vers Firebase"""
        if not self.firebase_initialized:
            return {'success': False, 'error': 'Firebase non initialisé'}
        
        conn = sqlite3.connect(self.local_db_path)
        cursor = conn.cursor()
        
        # Récupère enregistrements non synchronisés
        if record_id:
            cursor.execute("SELECT * FROM mesures_morpho WHERE id = ?", (record_id,))
        else:
            cursor.execute("SELECT * FROM mesures_morpho WHERE synced = 0")
        
        records = cursor.fetchall()
        results = {'success': 0, 'failed': 0, 'errors': []}
        
        for record in records:
            try:
                id_, animal_id, date_mesure, data_json, synced, cloud_id = record
                
                # Upload Firestore
                data = json.loads(data_json)
                doc_ref = self.db.collection('mesures_ovins').document()
                doc_ref.set({
                    'local_id': id_,
                    'animal_id': animal_id,
                    'date_mesure': date_mesure,
                    'data': data,
                    'synced_at': firestore.SERVER_TIMESTAMP
                })
                
                # Marque comme synchronisé
                cursor.execute(
                    "UPDATE mesures_morpho SET synced = 1, cloud_id = ? WHERE id = ?",
                    (doc_ref.id, id_)
                )
                results['success'] += 1
                
            except Exception as e:
                results['failed'] += 1
                results['errors'].append(str(e))
        
        conn.commit()
        conn.close()
        return results
    
    def sync_from_firebase(self, user_id: str = None) -> List[Dict]:
        """Récupère données cloud vers local"""
        if not self.firebase_initialized:
            return []
        
        try:
            query = self.db.collection('mesures_ovins')
            if user_id:
                query = query.where('user_id', '==', user_id)
            
            docs = query.order_by('date_mesure', direction=firestore.Query.DESCENDING).limit(100).stream()
            
            cloud_data = []
            for doc in docs:
                data = doc.to_dict()
                data['firestore_id'] = doc.id
                cloud_data.append(data)
            
            return cloud_data
            
        except Exception as e:
            st.error(f"Erreur récupération cloud: {e}")
            return []
    
    def backup_image_cloud(self, image_path: str, animal_id: str) -> str:
        """Upload image vers cloud storage"""
        if not self.firebase_initialized:
            return None
        
        try:
            blob_name = f"images/{animal_id}/{datetime.now().strftime('%Y%m%d_%H%M%S')}.jpg"
            blob = self.bucket.blob(blob_name)
            blob.upload_from_filename(image_path)
            blob.make_public()
            return blob.public_url
            
        except Exception as e:
            st.error(f"Erreur upload image: {e}")
            return None
    
    def render_cloud_interface(self):
        """Interface Streamlit synchronisation"""
        st.markdown("## ☁️ Synchronisation Cloud")
        
        # État connexion
        col_status1, col_status2 = st.columns(2)
        with col_status1:
            st.metric("Firebase", "✅ Connecté" if self.firebase_initialized else "❌ Déconnecté")
        with col_status2:
            st.metric("AWS S3", "✅ Connecté" if self.aws_client else "❌ Déconnecté")
        
        # Configuration
        with st.expander("🔧 Configuration Cloud"):
            cloud_provider = st.selectbox("Fournisseur", ["Firebase (Google)", "AWS S3"])
            
            if cloud_provider == "Firebase (Google)":
                st.info("Configuration via Streamlit Secrets (recommandé) ou fichier JSON")
                
                # Option 1: Secrets
                if st.checkbox("Utiliser secrets Streamlit"):
                    try:
                        firebase_secrets = st.secrets["firebase"]
                        if st.button("Initialiser Firebase"):
                            success = self.init_firebase(dict(firebase_secrets))
                            if success:
                                st.success("Firebase connecté!")
                    except:
                        st.error("Secrets Firebase non configurés dans .streamlit/secrets.toml")
                
                # Option 2: Fichier upload
                uploaded_cred = st.file_uploader("Fichier credentials JSON", type=['json'])
                if uploaded_cred:
                    cred_dict = json.load(uploaded_cred)
                    if st.button("Initialiser avec fichier"):
                        success = self.init_firebase(cred_dict)
                        if success:
                            st.success("Firebase connecté!")
            
            else:  # AWS
                aws_key = st.text_input("AWS Access Key", type="password")
                aws_secret = st.text_input("AWS Secret Key", type="password")
                aws_region = st.selectbox("Région", ["eu-west-1", "eu-west-3", "us-east-1"])
                
                if st.button("Connecter AWS"):
                    success = self.init_aws(aws_key, aws_secret, aws_region)
                    if success:
                        st.success("AWS S3 connecté!")
        
        # Synchronisation
        if self.firebase_initialized or self.aws_client:
            st.markdown("---")
            st.subheader("🔄 Gestion Synchronisation")
            
            col_sync1, col_sync2 = st.columns(2)
            
            with col_sync1:
                st.markdown("**Upload vers Cloud**")
                
                # Compteurs locaux
                conn = sqlite3.connect(self.local_db_path)
                cursor = conn.cursor()
                cursor.execute("SELECT COUNT(*) FROM mesures_morpho WHERE synced = 0")
                nb_unsynced = cursor.fetchone()[0]
                conn.close()
                
                st.metric("Enregistrements locaux non synchronisés", nb_unsynced)
                
                if nb_unsynced > 0 and st.button("☁️ Synchroniser maintenant", type="primary"):
                    with st.spinner("Synchronisation..."):
                        result = self.sync_to_firebase()
                        if result['success'] > 0:
                            st.success(f"✅ {result['success']} enregistrements synchronisés")
                        if result['failed'] > 0:
                            st.error(f"❌ {result['failed']} échecs")
            
            with col_sync2:
                st.markdown("**Download depuis Cloud**")
                
                if st.button("📥 Récupérer données cloud"):
                    with st.spinner("Téléchargement..."):
                        cloud_data = self.sync_from_firebase()
                        st.success(f"{len(cloud_data)} enregistrements récupérés")
                        
                        if cloud_data:
                            st.dataframe(pd.DataFrame([d['data']['mesure'] for d in cloud_data[:5]]))
        
        # Gestion offline
        st.markdown("---")
        st.subheader("💾 Gestion Offline")
        
        # Liste données locales
        conn = sqlite3.connect(self.local_db_path)
        df_local = pd.read_sql_query(
            "SELECT id, animal_id, date_mesure, synced FROM mesures_morpho ORDER BY date_mesure DESC", 
            conn
        )
        conn.close()
        
        if not df_local.empty:
            st.dataframe(df_local)
            
            # Export backup local
            if st.button("📦 Exporter backup SQLite"):
                with open(self.local_db_path, "rb") as f:
                    st.download_button("Télécharger .db", f.read(), 
                                     f"backup_ovin_{datetime.now().strftime('%Y%m%d')}.db",
                                     mime="application/x-sqlite3")
        else:
            st.info("Aucune donnée locale. Les mesures saisies seront stockées ici automatiquement.")

# ============================================================================
# INTÉGRATION DANS INTERFACE PRINCIPALE
# ============================================================================

def main():
    """Application principale complète"""
    st.set_page_config(
        page_title="Expert Ovin DZ Pro - Ultimate",
        page_icon="🐑",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Initialisation gestionnaires
    if 'cloud_manager' not in st.session_state:
        st.session_state.cloud_manager = CloudSyncManager()
    if 'pdf_manager' not in st.session_state:
        st.session_state.pdf_manager = PDFTemplateManager()
    if 'r_manager' not in st.session_state:
        st.session_state.r_manager = RExportManager()
    
    # Sidebar navigation
    st.sidebar.title("🐑 Expert Ovin DZ Pro")
    st.sidebar.caption("Version Ultimate - R + PDF + Cloud")
    
    module = st.sidebar.radio(
        "Navigation",
        [
            "🏠 Accueil",
            "📱 Analyse Image (Auto)",
            "📏 Saisie Manuelle Ruban",
            "📄 Templates PDF Terrain",      # NOUVEAU
            "📊 Statistiques R",              # NOUVEAU (avec R)
            "☁️ Synchronisation Cloud",       # NOUVEAU
            "🧬 Génétique & AlphaMissense",
            "🥛 Calcul IAL Complet",
            "⚙️ Configuration"
        ]
    )
    
    # Routage modules
    if module == "🏠 Accueil":
        st.title("🐑 Expert Ovin DZ Pro - Ultimate")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Données locales", len(st.session_state.get('mesures_manuelles', [])))
        with col2:
            st.metric("Module R", "✅ Actif" if RPY2_AVAILABLE else "⚠️ CSV only")
        with col3:
            st.metric("Cloud", "✅ Configurable")
        
        st.markdown("""
        ### Modules disponibles:
        
        1. **📱 Analyse Image** - Vision par ordinateur (OpenCV)
        2. **📏 Saisie Manuelle** - Protocole ruban standardisé
        3. **📄 PDF Terrain** - Fiches et guides imprimables
        4. **📊 Stats R** - Analyses statistiques avancées (rpy2)
        5. **☁️ Cloud** - Sync Firebase/AWS avec offline fallback
        6. **🧬 Génétique** - APIs NCBI/Ensembl/AlphaMissense
        7. **🥛 IAL** - Indice combiné Image/ADN/Pedigree
        """)
        
        # Derniers enregistrements
        if 'mesures_manuelles' in st.session_state and st.session_state.mesures_manuelles:
            st.subheader("Derniers enregistrements")
            dernieres = st.session_state.mesures_manuelles[-3:]
            for item in dernieres:
                m = item['mesure']
                st.write(f"🐑 {m.animal_id} ({m.race}) - {m.date_mesure.strftime('%d/%m/%Y')}")
    
    elif module == "📱 Analyse Image (Auto)":
        st.info("Module analyse image - Utiliser calibration automatique")
        # Appel classe ImageAnalyzerStreamlit existante
        
    elif module == "📏 Saisie Manuelle Ruban":
        from dataclasses import dataclass
        
        @dataclass
        class MesuresManuelles:
            animal_id: str
            date_mesure: datetime
            operateur: str
            race: str
            age_mois: int
            numero_lactation: int
            hauteur_garrot: float
            longueur_corps: float
            largeur_bassin: float
            tour_poitrine: float
            profondeur_poitrine: float
            largeur_epaules: float
            longueur_mamelle: float
            largeur_mamelle: float
            profondeur_mamelle: float
            circonference_tetine_gauche: float
            circonference_tetine_droite: float
            ecart_tetines: float
            hauteur_attache_mamelle: float
            angle_attache_mamelle: float
            angle_tetines: float
            score_attache_mamelle: int
            score_profondeur_mamelle: int
            score_symetrie: int
            score_etat_corps: int
            production_lait_jour: Optional[float] = None
            pourcentage_mg: Optional[float] = None
            pourcentage_mp: Optional[float] = None
            
            def to_dict(self):
                return {
                    'animal_id': self.animal_id,
                    'date_mesure': self.date_mesure.isoformat() if isinstance(self.date_mesure, datetime) else self.date_mesure,
                    'operateur': self.operateur,
                    'race': self.race,
                    'age_mois': self.age_mois,
                    'numero_lactation': self.numero_lactation,
                    'hauteur_garrot': self.hauteur_garrot,
                    'longueur_corps': self.longueur_corps,
                    'largeur_bassin': self.largeur_bassin,
                    'tour_poitrine': self.tour_poitrine,
                    'profondeur_poitrine': self.profondeur_poitrine,
                    'largeur_epaules': self.largeur_epaules,
                    'longueur_mamelle': self.longueur_mamelle,
                    'largeur_mamelle': self.largeur_mamelle,
                    'profondeur_mamelle': self.profondeur_mamelle,
                    'circonference_tetine_gauche': self.circonference_tetine_gauche,
                    'circonference_tetine_droite': self.circonference_tetine_droite,
                    'ecart_tetines': self.ecart_tetines,
                    'hauteur_attache_mamelle': self.hauteur_attache_mamelle,
                    'angle_attache_mamelle': self.angle_attache_mamelle,
                    'angle_tetines': self.angle_tetines,
                    'score_attache_mamelle': self.score_attache_mamelle,
                    'score_profondeur_mamelle': self.score_profondeur_mamelle,
                    'score_symetrie': self.score_symetrie,
                    'score_etat_corps': self.score_etat_corps,
                    'production_lait_jour': self.production_lait_jour,
                    'pourcentage_mg': self.pourcentage_mg,
                    'pourcentage_mp': self.pourcentage_mp
                }
            
            def calculer_indices(self):
                indices = {}
                if self.hauteur_garrot > 0:
                    indices['indice_conformation'] = round((self.longueur_corps / self.hauteur_garrot) * 100, 2)
                if self.tour_poitrine > 0:
                    poids_estime = (self.tour_poitrine ** 2) * self.longueur_corps / 10800
                    indices['poids_estime_kg'] = round(poids_estime, 2)
                    indices['indice_compacite'] = round((self.tour_poitrine / self.longueur_corps) * 100, 2)
                volume_mamelle = (4/3) * 3.14159 * (self.longueur_mamelle/2) * (self.largeur_mamelle/2) * (self.profondeur_mamelle/2)
                indices['volume_mamelle_cm3'] = round(volume_mamelle, 2)
                surface_mamelle = 3.14159 * (self.longueur_mamelle/2) * (self.largeur_mamelle/2)
                indices['surface_mamelle_cm2'] = round(surface_mamelle, 2)
                score_global = (
                    self.score_attache_mamelle * 0.3 +
                    self.score_profondeur_mamelle * 0.25 +
                    self.score_symetrie * 0.25 +
                    (10 - abs(self.score_etat_corps - 5)) * 0.2
                )
                indices['score_morphologique_global'] = round(score_global, 2)
                if self.circonference_tetine_droite > 0:
                    ratio = self.circonference_tetine_gauche / self.circonference_tetine_droite
                    indices['ratio_symetrie_tetines'] = round(ratio, 3)
                    indices['differentiel_tetines_cm'] = round(abs(self.circonference_tetine_gauche - self.circonference_tetine_droite), 2)
                return indices
        
        class SaisieManuelleRuban:
            def __init__(self):
                self.session_state = st.session_state
                if 'mesures_manuelles' not in self.session_state:
                    self.session_state.mesures_manuelles = []
            
            def render_interface(self):
                st.markdown("## 📏 Saisie Manuelle des Mesures (Ruban Métrique)")
                st.info("""
                **Protocole de mesure standardisé (ISO 7481 - Ovin)**
                
                1. Animal en station debout, sur sol plat
                2. Ruban métrique souple, contact sans compression
                3. Mesures effectuées par opérateur formé
                4. Noter la date et l'heure de la mesure
                """)
                
                with st.form("saisie_mesures_manuelles"):
                    st.subheader("🐑 Identification")
                    
                    col_id1, col_id2, col_id3 = st.columns(3)
                    with col_id1:
                        animal_id = st.text_input("ID Animal", f"MOUTON_{datetime.now().strftime('%Y%m%d_%H%M')}")
                        race = st.selectbox("Race", ["Lacaune", "Manech", "Basco-Béarnaise", "Corsican", "Awassi", "Dorper", "Autre"])
                    with col_id2:
                        date_mesure = st.date_input("Date mesure", datetime.now())
                        age_mois = st.number_input("Âge (mois)", 8, 180, 24)
                    with col_id3:
                        operateur = st.text_input("Opérateur", "Technicien")
                        numero_lactation = st.number_input("N° lactation", 0, 10, 1)
                    
                    st.markdown("---")
                    st.subheader("📐 Mensurations Corporelles (cm)")
                    
                    col_corp1, col_corp2, col_corp3 = st.columns(3)
                    with col_corp1:
                        hauteur_garrot = st.number_input("Hauteur au garrot", 50.0, 90.0, 68.0, 0.5, help="Point le plus haut du dos, verticale jusqu'au sol")
                        longueur_corps = st.number_input("Longueur du corps", 60.0, 100.0, 78.0, 0.5, help="Point d'épaule à pointe fessière")
                        largeur_epaules = st.number_input("Largeur des épaules", 15.0, 30.0, 22.0, 0.5)
                    
                    with col_corp2:
                        largeur_bassin = st.number_input("Largeur du bassin", 15.0, 30.0, 21.0, 0.5, help="Distance entre tubers coxaux")
                        tour_poitrine = st.number_input("Tour de poitrine", 70.0, 120.0, 92.0, 0.5, help="Arrière des épaules autour du poitrail")
                        profondeur_poitrine = st.number_input("Profondeur poitrine", 25.0, 45.0, 34.0, 0.5, help="Garrot au sternum derrière les épaules")
                    
                    with col_corp3:
                        st.markdown("**Références race Lacaune**")
                        st.caption(f"Garrot optimal: 68 cm")
                        st.caption(f"Longueur optimale: 78 cm")
                        st.progress(min(1.0, 68.0/90))
                        st.caption("Hauteur (normalisée)")
                        st.progress(min(1.0, 78.0/100))
                        st.caption("Longueur (normalisée)")
                    
                    st.markdown("---")
                    st.subheader("🥛 Mensurations Mamelle (cm)")
                    
                    col_mam1, col_mam2, col_mam3 = st.columns(3)
                    with col_mam1:
                        longueur_mamelle = st.number_input("Longueur mamelle", 10.0, 30.0, 18.0, 0.5, help="Jonction corps à pointe arrière")
                        largeur_mamelle = st.number_input("Largeur mamelle", 8.0, 25.0, 15.0, 0.5, help="Maximum au niveau des tétines")
                        profondeur_mamelle = st.number_input("Profondeur mamelle", 8.0, 25.0, 16.0, 0.5, help="Attache à pointe inférieure")
                    
                    with col_mam2:
                        circonf_gauche = st.number_input("Circonf. tétine G", 2.0, 5.0, 3.2, 0.1, help="Base de la tétine")
                        circonf_droite = st.number_input("Circonf. tétine D", 2.0, 5.0, 3.2, 0.1)
                        ecart_tetines = st.number_input("Écart tétines", 5.0, 20.0, 11.0, 0.5, help="Distance entre centres des tétines")
                    
                    with col_mam3:
                        hauteur_attache = st.number_input("Hauteur attache", 20.0, 50.0, 35.0, 0.5, help="Du sol à l'attache antérieure")
                        angle_attache = st.number_input("Angle attache (°)", 30.0, 90.0, 60.0, 1.0, help="Angle avec la paroi abdominale")
                        angle_tetines = st.number_input("Angle tétines (°)", 0.0, 45.0, 15.0, 1.0, help="Inclinaison vers l'extérieur")
                    
                    st.markdown("---")
                    st.subheader("⭐ Scores Subjectifs (échelle 1-9)")
                    
                    col_score1, col_score2, col_score3, col_score4 = st.columns(4)
                    with col_score1:
                        score_attache = st.slider("Attache mamelle", 1, 9, 7, help="1=Très lâche, 9=Très serrée")
                    with col_score2:
                        score_profondeur = st.slider("Profondeur mamelle", 1, 9, 7, help="1=Superficielle, 9=Très profonde")
                    with col_score3:
                        score_symetrie = st.slider("Symétrie", 1, 9, 8, help="1=Asymétrique, 9=Parfaite")
                    with col_score4:
                        score_ec = st.slider("État corporel", 1, 9, 5, help="1=Émacié, 5=Optimal, 9=Obèse")
                    
                    st.markdown("---")
                    st.subheader("🥛 Production Laitière (optionnel)")
                    
                    col_prod1, col_prod2, col_prod3 = st.columns(3)
                    with col_prod1:
                        prod_lait = st.number_input("Production/jour (L)", 0.0, 5.0, 0.0, 0.1) or None
                    with col_prod2:
                        pc_mg = st.number_input("% Matière Grasse", 0.0, 15.0, 0.0, 0.1) or None
                    with col_prod3:
                        pc_mp = st.number_input("% Matière Protéique", 0.0, 10.0, 0.0, 0.1) or None
                    
                    submitted = st.form_submit_button("💾 Enregistrer les mesures", type="primary")
                    
                    if submitted:
                        mesure = MesuresManuelles(
                            animal_id=animal_id,
                            date_mesure=datetime.combine(date_mesure, datetime.min.time()),
                            operateur=operateur,
                            race=race,
                            age_mois=age_mois,
                            numero_lactation=numero_lactation,
                            hauteur_garrot=hauteur_garrot,
                            longueur_corps=longueur_corps,
                            largeur_bassin=largeur_bassin,
                            tour_poitrine=tour_poitrine,
                            profondeur_poitrine=profondeur_poitrine,
                            largeur_epaules=largeur_epaules,
                            longueur_mamelle=longueur_mamelle,
                            largeur_mamelle=largeur_mamelle,
                            profondeur_mamelle=profondeur_mamelle,
                            circonference_tetine_gauche=circonf_gauche,
                            circonference_tetine_droite=circonf_droite,
                            ecart_tetines=ecart_tetines,
                            hauteur_attache_mamelle=hauteur_attache,
                            angle_attache_mamelle=angle_attache,
                            angle_tetines=angle_tetines,
                            score_attache_mamelle=score_attache,
                            score_profondeur_mamelle=score_profondeur,
                            score_symetrie=score_symetrie,
                            score_etat_corps=score_ec,
                            production_lait_jour=prod_lait if prod_lait > 0 else None,
                            pourcentage_mg=pc_mg if pc_mg > 0 else None,
                            pourcentage_mp=pc_mp if pc_mp > 0 else None
                        )
                        
                        indices = mesure.calculer_indices()
                        
                        # Sauvegarde session + cloud local
                        record = {'mesure': mesure, 'indices': indices}
                        self.session_state.mesures_manuelles.append(record)
                        
                        # Sauvegarde SQLite (via cloud manager)
                        st.session_state.cloud_manager.save_measurement_local(mesure, indices)
                        
                        st.success(f"✅ Mesures enregistrées pour {animal_id} (ID local: {st.session_state.cloud_manager.local_db_path})")
                        
                        col_res1, col_res2 = st.columns(2)
                        with col_res1:
                            st.metric("Poids estimé", f"{indices.get('poids_estime_kg', 0)} kg")
                            st.metric("Indice conformation", f"{indices.get('indice_conformation', 0):.1f}")
                        with col_res2:
                            st.metric("Volume mamelle", f"{indices.get('volume_mamelle_cm3', 0)} cm³")
                            st.metric("Score global", f"{indices.get('score_morphologique_global', 0)}/9")
                        
                        if indices.get('indice_conformation', 0) < 100:
                            st.warning("⚠️ Indice de conformation bas - Vérifier mensurations")
                        if indices.get('differentiel_tetines_cm', 0) > 0.5:
                            st.warning("⚠️ Différence significative entre tétines - Asymétrie à noter")
        
        saisie = SaisieManuelleRuban()
        saisie.render_interface()
        
    elif module == "📄 Templates PDF Terrain":
        st.session_state.pdf_manager.render_pdf_interface()
        
    elif module == "📊 Statistiques R":
        # Récupération données
        mesures = st.session_state.get('mesures_manuelles', [])
        
        if mesures:
            # Conversion en DataFrame pour R
            df_data = []
            for item in mesures:
                m = item['mesure']
                row = m.to_dict()
                row.update(item['indices'])
                df_data.append(row)
            
            df = pd.DataFrame(df_data)
            
            # Interface R
            st.session_state.r_manager.render_r_export_interface(df)
            
            # Interface stats Python (fallback)
            with st.expander("Stats Python (fallback si R indisponible)"):
                st.dataframe(df.describe())
        else:
            st.warning("Aucune donnée disponible. Veuillez d'abord saisir des mesures.")
        
    elif module == "☁️ Synchronisation Cloud":
        st.session_state.cloud_manager.render_cloud_interface()
        
    elif module == "🧬 Génétique & AlphaMissense":
        st.info("Module génétique - NCBI/Ensembl/AlphaMissense")
        
    elif module == "🥛 Calcul IAL Complet":
        st.info("Module IAL - Intégration toutes sources")
        
    elif module == "⚙️ Configuration":
        st.title("⚙️ Configuration Système")
        
        st.subheader("Dépendances détectées")
        deps = {
            'rpy2 (Export R)': RPY2_AVAILABLE,
            'reportlab (PDF)': REPORTLAB_AVAILABLE,
            'firebase-admin (Cloud)': FIREBASE_AVAILABLE,
            'boto3 (AWS)': AWS_AVAILABLE,
            'opencv-python (Vision)': OPENCV_AVAILABLE,
            'statsmodels (Stats)': STATSMODELS_AVAILABLE,
            'scikit-learn (ML)': SKLEARN_AVAILABLE
        }
        
        for lib, status in deps.items():
            st.write(f"{'✅' if status else '❌'} {lib}")
        
        st.subheader("Configuration Stockage Local")
        st.write(f"Base SQLite: `{st.session_state.cloud_manager.local_db_path}`")

if __name__ == "__main__":
    main()
