"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D et Génétique
Base de données simulée de races ovines algériennes
"""

# ========== IMPORTS ==========
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
import random
import math
import io
import base64
from PIL import Image, ImageDraw

# ========== CONFIGURATION STREAMLIT ==========
st.set_page_config(
    page_title="Ovin Manager Pro - Races Algériennes",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ========== CSS PERSONNALISÉ ==========
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
</style>
""", unsafe_allow_html=True)

# ========== STANDARDS DES RACES ALGÉRIENNES ==========
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
        'prolificite': (1.2, 1.8)
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
        'prolificite': (1.1, 1.5)
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
        'prolificite': (1.3, 1.7)
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
        'prolificite': (1.0, 1.4)
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
        'prolificite': (1.2, 1.8)
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
        'prolificite': (1.0, 1.6)
    }
}

# ========== FONCTIONS STATISTIQUES (sans scipy) ==========
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

# ========== BASE DE DONNÉES ==========
def creer_base_races():
    """Crée une base de données avec races algériennes"""
    conn = sqlite3.connect('ovin_algerien.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table brebis
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE NOT NULL,
            nom TEXT,
            race TEXT CHECK(race IN ('HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU')),
            sous_race TEXT,
            sexe TEXT CHECK(sexe IN ('F', 'M')),
            date_naissance DATE,
            age_mois INTEGER,
            poids FLOAT,
            score_condition INTEGER CHECK(score_condition BETWEEN 1 AND 5),
            couleur_robe TEXT,
            intensite_couleur INTEGER CHECK(intensite_couleur BETWEEN 1 AND 10),
            cornes BOOLEAN,
            taille_cornes_cm FLOAT,
            forme_cornes TEXT,
            type_laine TEXT,
            qualite_laine INTEGER CHECK(qualite_laine BETWEEN 1 AND 10),
            longueur_corps_cm FLOAT,
            hauteur_garrot_cm FLOAT,
            largeur_bassin_cm FLOAT,
            tour_poitrine_cm FLOAT,
            circonference_tete_cm FLOAT,
            longueur_oreille_cm FLOAT,
            temperement TEXT CHECK(temperement IN ('calme', 'nervieux', 'intermediaire')),
            aptitude TEXT CHECK(aptitude IN ('lait', 'viande', 'mixte', 'laine')),
            notes TEXT,
            mere_id TEXT,
            pere_id TEXT,
            coefficient_consanguinite FLOAT DEFAULT 0.0,
            statut TEXT DEFAULT 'active',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Tables supplémentaires
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
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
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
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
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
            date_analyse DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Peupler la base si vide
    cursor.execute("SELECT COUNT(*) FROM brebis")
    if cursor.fetchone()[0] == 0:
        peupler_base_races(cursor, conn)
    
    conn.commit()
    return conn

def peupler_base_races(cursor, conn):
    """Peuple la base avec des races algériennes"""
    races = ['HAMRA', 'OUDA', 'SIDAHOU', 'BERBERE', 'CROISE', 'INCONNU']
    brebis_data = []
    
    for i in range(1, 51):
        race = random.choice(races)
        sexe = random.choices(['F', 'M'], weights=[0.7, 0.3])[0]
        
        # Générer identifiant
        identifiant = f"{race[:3]}-{sexe}-2023-{i:03d}"
        nom = f"Brebis {i}"
        age_mois = random.randint(12, 84)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        
        # Poids selon race et sexe
        poids_min, poids_max = STANDARDS_RACES[race]['poids_adulte'][sexe.lower()]
        poids = random.uniform(poids_min, poids_max)
        
        brebis_data.append((
            identifiant, nom, race, '', sexe, date_naissance.isoformat(), 
            age_mois, poids, random.randint(2, 4), 'Couleur', 
            random.randint(5, 10), random.choice([True, False]), 
            random.uniform(0, 60), '', 'type', random.randint(3, 9),
            random.uniform(*STANDARDS_RACES[race]['mensurations']['longueur_cm']),
            random.uniform(*STANDARDS_RACES[race]['mensurations']['hauteur_cm']),
            random.uniform(*STANDARDS_RACES[race]['mensurations']['largeur_bassin_cm']),
            random.uniform(*STANDARDS_RACES[race]['mensurations']['tour_poitrine_cm']),
            random.uniform(45, 65), random.uniform(12, 18),
            random.choice(['calme', 'nervieux', 'intermediaire']),
            random.choice(['lait', 'viande', 'mixte', 'laine']),
            f"Brebis {race} - Élevage algérien",
            None, None, random.uniform(0.0, 0.15), 'active'
        ))
    
    cursor.executemany('''
        INSERT INTO brebis (
            identifiant, nom, race, sous_race, sexe, date_naissance, age_mois, 
            poids, score_condition, couleur_robe, intensite_couleur, cornes, 
            taille_cornes_cm, forme_cornes, type_laine, qualite_laine,
            longueur_corps_cm, hauteur_garrot_cm, largeur_bassin_cm, 
            tour_poitrine_cm, circonference_tete_cm, longueur_oreille_cm,
            temperement, aptitude, notes, mere_id, pere_id, 
            coefficient_consanguinite, statut
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', brebis_data)
    
    conn.commit()

# Initialiser la base
conn = creer_base_races()

# ========== PAGE CRITÈRES DE SÉLECTION ==========
def page_criteres():
    """Page des critères de sélection morphologiques et phénotypiques"""
    st.markdown('<h2 class="section-header">🎯 CRITÈRES DE SÉLECTION - MAMMELLES</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["📏 MAMMELLES", "🏋️ MORPHOLOGIE", "🧬 PHÉNOTYPE", "📊 SCORING"])
    
    with tab1:
        st.markdown("### 📏 CRITÈRES MAMMAIRES - PRODUCTION LAITIÈRE")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            <div class='metric-card'>
                <h4>🎯 PRINCIPAUX CRITÈRES</h4>
                <p><strong>1. Volume mammaire</strong></p>
                <p><strong>2. Insertion des trayons</strong></p>
                <p><strong>3. Symétrie</strong></p>
                <p><strong>4. Longueur des trayons</strong></p>
                <p><strong>5. Orientation</strong></p>
            </div>
            """, unsafe_allow_html=True)
            
            # Formulaire d'évaluation des mamelles
            with st.form("evaluation_mamelles"):
                st.markdown("#### 📝 ÉVALUATION MAMMAIRE")
                
                volume = st.slider("Volume mammaire (1-5)", 1, 5, 3, 
                                  help="1: Très petit, 5: Très développé")
                symetrie = st.slider("Symétrie (1-5)", 1, 5, 3,
                                    help="1: Asymétrique, 5: Parfaitement symétrique")
                insertion = st.slider("Insertion des trayons (1-5)", 1, 5, 3,
                                     help="1: Très écartés, 5: Bien insérés")
                longueur_trayons = st.slider("Longueur des trayons (cm)", 2.0, 8.0, 4.5, 0.1)
                
                if st.form_submit_button("💾 Évaluer", type="primary"):
                    score_total = (volume + symetrie + insertion) / 3
                    st.success(f"✅ Score mammelle: {score_total:.1f}/5")
        
        with col2:
            # Schéma explicatif
            st.markdown("""
            <div style='text-align: center; padding: 20px; background: #f8f9fa; border-radius: 10px;'>
                <h4>🏷️ CLASSIFICATION MAMMAIRE</h4>
                <p><strong>Type A:</strong> Mamelle bien développée, trayons parallèles</p>
                <p><strong>Type B:</strong> Mamelle moyenne, légère asymétrie</p>
                <p><strong>Type C:</strong> Mamelle faible, trayons divergents</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Critères spécifiques races
            st.markdown("#### 🏷️ STANDARDS PAR RACE")
            
            criteres_races = {
                'HAMRA': {
                    'volume': 'Moyen à élevé',
                    'trayons': '3-5 cm, bien orientés',
                    'capacite': 'Bonne adaptabilité désert'
                },
                'OUDA': {
                    'volume': 'Grand',
                    'trayons': '4-6 cm, légèrement divergents',
                    'capacite': 'Excellente production'
                },
                'BERBERE': {
                    'volume': 'Petit à moyen',
                    'trayons': '2-4 cm, bien insérés',
                    'capacite': 'Adaptée montagne'
                }
            }
            
            race_selection = st.selectbox("Voir standards", list(criteres_races.keys()))
            if race_selection in criteres_races:
                st.info(f"""
                **{race_selection} - Critères mammaires:**
                - Volume: {criteres_races[race_selection]['volume']}
                - Trayons: {criteres_races[race_selection]['trayons']}
                - Capacité: {criteres_races[race_selection]['capacite']}
                """)
    
    with tab2:
        st.markdown("### 🏋️ CRITÈRES MORPHOLOGIQUES GÉNÉRAUX")
        
        col_morph1, col_morph2 = st.columns(2)
        
        with col_morph1:
            st.markdown("""
            <div class='metric-card'>
                <h4>📏 MENSURATIONS</h4>
                <p><strong>• Longueur corps:</strong> 90-130 cm</p>
                <p><strong>• Hauteur garrot:</strong> 60-90 cm</p>
                <p><strong>• Largeur bassin:</strong> 35-55 cm</p>
                <p><strong>• Tour poitrine:</strong> 95-135 cm</p>
                <p><strong>• Indice corporel:</strong> 2.5-3.5</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Indices de conformation
            st.markdown("#### 📊 INDICES DE CONFORMATION")
            
            indices = {
                'Indice de format': 'Longueur / Hauteur × 100',
                'Indice thoracique': 'Tour poitrine / Hauteur × 100',
                'Indice de bassin': 'Largeur bassin / Hauteur × 100',
                'Indice laitier': 'Tour poitrine / Poids × 100'
            }
            
            for nom, formule in indices.items():
                with st.expander(f"{nom}"):
                    st.write(f"**Formule:** {formule}")
                    st.write("**Interprétation:**")
                    if "laitier" in nom.lower():
                        st.write("> 65: Type laitier, < 55: Type viande")
        
        with col_morph2:
            st.markdown("""
            <div class='metric-card'>
                <h4>🎨 CARACTÈRES QUALITATIFS</h4>
                <p><strong>• Robe:</strong> Couleur, intensité, uniformité</p>
                <p><strong>• Cornes:</strong> Présence, forme, taille</p>
                <p><strong>• Laine:</strongType, finesse, densité</p>
                <p><strong>• Tête:</strong> Forme, profil, oreilles</p>
                <p><strong>• Membres:</strong> Solides, bien d'aplomb</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Système de scoring
            st.markdown("#### 🏆 SYSTÈME DE SCORING")
            
            caracteres = st.multiselect(
                "Caractères à évaluer:",
                ['Conformation générale', 'Développement musculaire', 'Qualité squelette', 
                 'Membres', 'Robe', 'Tempérament', 'Aptitude laitière'],
                default=['Conformation générale', 'Aptitude laitière']
            )
            
            if st.button("📊 Calculer score total"):
                score = len(caracteres) * 3  # Score simulé
                st.metric("Score conformation", f"{score}/21", "Bon")
    
    with tab3:
        st.markdown("### 🧬 CRITÈRES PHÉNOTYPIQUES")
        
        # Phénotypes liés à la production
        st.markdown("#### 🥛 PHÉNOTYPES LAITIERS")
        
        phenotypes = pd.DataFrame({
            'Caractère': ['Développement mammaire', 'Veines mammaires', 'Tempérament', 
                         'Appétit', 'Rapidité de traite', 'Résistance mammite'],
            'Héritabilité': [0.35, 0.25, 0.15, 0.20, 0.30, 0.10],
            'Impact production': ['Élevé', 'Moyen', 'Faible', 'Moyen', 'Élevé', 'Élevé']
        })
        
        st.dataframe(phenotypes.style.background_gradient(subset=['Héritabilité']))
        
        # Corrélations phénotypiques
        st.markdown("#### 🔗 CORRÉLATIONS PHÉNOTYPIQUES")
        
        fig = go.Figure(data=go.Heatmap(
            z=[[1, 0.6, 0.3, 0.4],
               [0.6, 1, 0.2, 0.5],
               [0.3, 0.2, 1, 0.1],
               [0.4, 0.5, 0.1, 1]],
            x=['Production lait', 'Taux MG', 'Taux protéine', 'Cellules'],
            y=['Production lait', 'Taux MG', 'Taux protéine', 'Cellules'],
            colorscale='RdBu',
            text=[[1, 0.6, 0.3, 0.4],
                  [0.6, 1, 0.2, 0.5],
                  [0.3, 0.2, 1, 0.1],
                  [0.4, 0.5, 0.1, 1]],
            texttemplate='%{text:.2f}',
            textfont={"size": 10}
        ))
        
        fig.update_layout(title="Corrélations entre caractères laitiers")
        st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        st.markdown("### 📊 SCORING INTÉGRÉ")
        
        # Formulaire complet de scoring
        with st.form("scoring_complet"):
            st.markdown("#### 🎯 ÉVALUATION COMPLÈTE")
            
            col_score1, col_score2, col_score3 = st.columns(3)
            
            with col_score1:
                st.markdown("**📏 MORPHOLOGIE**")
                conformation = st.slider("Conformation (1-10)", 1, 10, 7)
                developpement = st.slider("Développement (1-10)", 1, 10, 6)
            
            with col_score2:
                st.markdown("**🥛 MAMMELLES**")
                volume_m = st.slider("Volume (1-10)", 1, 10, 8)
                symetrie_m = st.slider("Symétrie (1-10)", 1, 10, 7)
                insertion_m = st.slider("Insertion (1-10)", 1, 10, 6)
            
            with col_score3:
                st.markdown("**🧬 GÉNÉTIQUE**")
                valeur_gen = st.slider("Valeur génétique (1-10)", 1, 10, 7)
                diversite = st.slider("Diversité (1-10)", 1, 10, 8)
            
            notes = st.text_area("Observations complémentaires")
            
            if st.form_submit_button("🎯 Calculer score final", type="primary"):
                # Calcul des scores
                score_morph = (conformation + developpement) / 2
                score_mamelle = (volume_m + symetrie_m + insertion_m) / 3
                score_gen = (valeur_gen + diversite) / 2
                score_final = (score_morph * 0.3 + score_mamelle * 0.4 + score_gen * 0.3)
                
                # Affichage résultats
                col_res1, col_res2, col_res3 = st.columns(3)
                
                with col_res1:
                    st.metric("Score morphologie", f"{score_morph:.1f}/10")
                
                with col_res2:
                    st.metric("Score mamelles", f"{score_mamelle:.1f}/10")
                
                with col_res3:
                    st.metric("Score génétique", f"{score_gen:.1f}/10")
                
                # Score final avec interprétation
                st.markdown("---")
                st.markdown(f"### 🏆 SCORE FINAL: **{score_final:.1f}/10**")
                
                if score_final >= 8:
                    st.success("🎖️ **EXCELLENT** - Animal d'élite pour la reproduction")
                elif score_final >= 6:
                    st.info("✅ **BON** - Animal de production satisfaisant")
                elif score_final >= 4:
                    st.warning("⚠️ **MOYEN** - À surveiller ou améliorer")
                else:
                    st.error("❌ **FAIBLE** - À réformer ou surveiller étroitement")
                
                # Recommandations
                if score_mamelle < 5:
                    st.warning("**Recommandation:** Améliorer la sélection sur les critères mammaires")
                if score_gen < 5:
                    st.warning("**Recommandation:** Introduire de nouveaux reproducteurs")

# ========== MODULE SCANNER 3D ==========
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
            'HAMRA': (139, 0, 0), 'OUDA': (255, 255, 255),
            'SIDAHOU': (50, 50, 50), 'BERBERE': (165, 42, 42),
            'CROISE': (160, 120, 80), 'INCONNU': (200, 200, 200)
        }
        
        race = brebis_info.get('race', 'INCONNU')
        corps_color = couleurs.get(race, (200, 200, 200))
        
        # Dessiner le corps
        draw.ellipse([100, 80, 300, 200], fill=corps_color, outline='black', width=2)
        draw.ellipse([280, 100, 350, 160], fill=corps_color, outline='black', width=2)
        
        # Pattes
        for x in [130, 170, 230, 270]:
            draw.rectangle([x, 200, x+20, 280], fill='black')
        
        # Informations
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {race}", fill='black')
        draw.text((10, 50), f"Poids: {brebis_info.get('poids', 0):.1f} kg", fill='black')
        
        return image

# ========== MODULE GÉNÉTIQUE ==========
class ModuleGenetique:
    """Module d'analyse génétique"""
    
    @staticmethod
    def generer_genotype(brebis_id, race):
        """Génère un génotype simulé"""
        genotypes = []
        
        for i in range(10):
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
        """Calcule la diversité génétique - VERSION CORRIGÉE"""
        if not genotypes:
            return {}
        
        # Convertir en DataFrame de manière sécurisée
        data = []
        for geno in genotypes:
            # S'assurer que nous avons assez d'éléments
            if len(geno) >= 8:  # Minimum pour les calculs de base
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
        heterozygotes = df[df['allele1'] != df['allele2']]
        ho = len(heterozygotes) / len(df) if len(df) > 0 else 0
        he = 1 - (df['freq_allelique']**2).mean() if 'freq_allelique' in df.columns else 0
        fis = 1 - (ho / he) if he > 0 else 0
        
        return {
            'heterozygosite_observee': round(ho, 4),
            'heterozygosite_attendue': round(he, 4),
            'fis': round(fis, 4),
            'nombre_snps': len(df['marqueur'].unique())
        }

# ========== PAGE ACCUEIL ==========
def page_accueil():
    """Page d'accueil avec vue d'ensemble"""
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - RACES ALGÉRIENNES</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique des races ovines algériennes**")
    
    # Métriques principales
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

# ========== PAGE SCANNER 3D ==========
def page_scanner_3d():
    """Page du scanner 3D avec saisie manuelle"""
    st.markdown('<h2 class="section-header">📐 SCANNER 3D & SAISIE MANUELLE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🎯 SCANNER 3D", "📝 SAISIE MANUELLE"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY nom")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_option = st.selectbox(
                "SÉLECTIONNEZ UNE BREBIS À SCANNER:",
                [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list]
            )
            
            if brebis_option:
                # Extraction sécurisée de l'ID
                try:
                    brebis_id = int(brebis_option.split('(')[1].split(')')[0].split('-')[-1])
                except:
                    st.error("Erreur dans le format de l'identifiant")
                    return
                
                # Photo simulée
                cursor.execute("SELECT * FROM brebis WHERE id = ?", (brebis_id,))
                brebis_info = cursor.fetchone()
                columns = [desc[0] for desc in cursor.description]
                brebis_dict = dict(zip(columns, brebis_info))
                
                photo = Scanner3D.generer_photo_simulee(brebis_dict)
                st.image(photo, caption=f"Photo simulée - {brebis_dict['nom']}", use_column_width=True)

# ========== PAGE GESTION ==========
def page_gestion():
    """Page de gestion du troupeau"""
    st.markdown('<h2 class="section-header">📊 GESTION DU TROUPEAU</h2>', unsafe_allow_html=True)
    
    cursor = conn.cursor()
    cursor.execute("""
        SELECT identifiant, nom, race, sexe, age_mois, poids, 
               score_condition, couleur_robe, statut
        FROM brebis
        ORDER BY race, identifiant
    """)
    
    brebis_data = cursor.fetchall()
    df = pd.DataFrame(brebis_data, columns=['ID', 'Nom', 'Race', 'Sexe', 'Âge', 'Poids', 'Score', 'Couleur', 'Statut'])
    
    st.dataframe(df, use_container_width=True, height=500)

# ========== PAGE PRODUCTION ==========
def page_production():
    """Page de suivi de production"""
    st.markdown('<h2 class="section-header">🥛 SUIVI DE PRODUCTION</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["📝 SAISIE", "📈 ANALYSE"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom, race FROM brebis WHERE sexe = 'F' AND statut = 'active'")
        femelles = cursor.fetchall()
        
        if femelles:
            with st.form("form_production"):
                brebis_sel = st.selectbox("Sélectionner une brebis", 
                                        [f"{f[1]} ({f[2]})" for f in femelles])
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5, 0.1)
                
                with col2:
                    mg = st.number_input("Matière grasse %", 0.0, 20.0, 7.2, 0.1)
                
                with col3:
                    proteine = st.number_input("Protéine %", 0.0, 20.0, 5.5, 0.1)
                
                date_mesure = st.date_input("Date", value=date.today())
                
                if st.form_submit_button("💾 Enregistrer", type="primary"):
                    st.success("✅ Production enregistrée!")

# ========== PAGE STATISTIQUES (RSTATS) ==========
def page_stats():
    """Page d'analyse statistique avancée"""
    st.markdown('<h2 class="section-header">📊 ANALYSE STATISTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    cursor = conn.cursor()
    cursor.execute("""
        SELECT race, sexe, age_mois, poids, score_condition, 
               longueur_corps_cm, hauteur_garrot_cm
        FROM brebis
        WHERE statut = 'active'
    """)
    
    data = cursor.fetchall()
    df = pd.DataFrame(data, columns=['Race', 'Sexe', 'Âge', 'Poids', 'Score', 'Longueur', 'Hauteur'])
    
    # Statistiques descriptives
    st.markdown("### 📈 STATISTIQUES DESCRIPTIVES")
    
    if not df.empty:
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        
        for col in numeric_cols:
            st.write(f"**{col}**")
            col1, col2, col3, col4 = st.columns(4)
            with col1: st.metric("Moyenne", f"{df[col].mean():.2f}")
            with col2: st.metric("Écart-type", f"{df[col].std():.2f}")
            with col3: st.metric("Minimum", f"{df[col].min():.2f}")
            with col4: st.metric("Maximum", f"{df[col].max():.2f}")

# ========== PAGE GÉNÉTIQUE ==========
def page_genetique():
    """Page d'analyse génétique avancée - VERSION CORRIGÉE"""
    st.markdown('<h2 class="section-header">🧬 ANALYSE GÉNÉTIQUE AVANCÉE</h2>', unsafe_allow_html=True)
    
    tab1, tab2 = st.tabs(["🧪 GÉNOTYPAGE", "🌳 DIVERSITÉ"])
    
    with tab1:
        st.markdown("### 🧪 GÉNOTYPAGE SNP")
        
        cursor = conn.cursor()
        cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY race")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_select = st.selectbox("Sélectionner une brebis", 
                                        [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list])
            
            if brebis_select and st.button("🧬 Générer génotype", type="primary"):
                try:
                    # Extraction sécurisée
                    brebis_id = int(brebis_select.split('(')[1].split(')')[0].split('-')[-1])
                    race = brebis_select.split('- ')[1]
                    
                    # Générer génotype
                    genotypes = ModuleGenetique.generer_genotype(brebis_id, race)
                    
                    # Insérer dans la base
                    for genotype in genotypes:
                        cursor.execute('''
                            INSERT OR REPLACE INTO genotypage 
                            (brebis_id, marqueur, chromosome, position, allele1, allele2, 
                             genotype, frequence_allelique, effet_additif, effet_dominant, 
                             r2, p_value, gene_associe, trait_associe, date_analyse)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', genotype)
                    
                    conn.commit()
                    st.success(f"✅ Génotype généré pour {brebis_select}")
                    
                except Exception as e:
                    st.error(f"Erreur: {str(e)}")
    
    with tab2:
        st.markdown("### 🌳 DIVERSITÉ GÉNÉTIQUE")
        
        cursor.execute("""
            SELECT g.brebis_id, b.race, b.identifiant,
                   g.marqueur, g.allele1, g.allele2, g.genotype,
                   g.frequence_allelique, g.trait_associe
            FROM genotypage g
            JOIN brebis b ON g.brebis_id = b.id
            LIMIT 50
        """)
        
        geno_data = cursor.fetchall()
        
        if geno_data:
            try:
                # UTILISER LA FONCTION CORRIGÉE
                diversite = ModuleGenetique.calculer_diversite_genetique(geno_data)
                
                if diversite:
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        st.metric("Hétérozygotie observée", f"{diversite['heterozygosite_observee']:.4f}")
                    with col2:
                        st.metric("Hétérozygotie attendue", f"{diversite['heterozygosite_attendue']:.4f}")
                    with col3:
                        st.metric("Fis", f"{diversite['fis']:.4f}")
                    with col4:
                        st.metric("SNPs analysés", diversite['nombre_snps'])
                else:
                    st.info("Impossible de calculer la diversité génétique")
                    
            except Exception as e:
                st.error(f"Erreur dans le calcul: {str(e)}")
                st.info("Génération de données de test pour démonstration...")
                
                # Données de test
                test_data = []
                for i in range(10):
                    test_data.append((
                        1, 'TEST', 'TEST001',
                        f'SNP{i:03d}', 'A', 'G', 'AG',
                        0.5, 'test'
                    ))
                
                diversite = ModuleGenetique.calculer_diversite_genetique(test_data)
                if diversite:
                    st.info(f"Diversité calculée sur données de test: H = {diversite['heterozygosite_observee']:.3f}")
        else:
            st.info("Aucune donnée de génotypage disponible. Générez d'abord des génotypes.")

# ========== BARRE LATÉRALE ==========
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
         "📐 SCANNER 3D", 
         "📊 GESTION", 
         "🥛 PRODUCTION",
         "🎯 CRITÈRES",
         "📊 RSTATS",
         "🧬 GÉNÉTIQUE"]
    )
    
    st.markdown("---")
    
    # Aide rapide sur les critères
    if page == "🎯 CRITÈRES":
        st.markdown("### 💡 AIDE RAPIDE")
        st.info("""
        **Critères mammaires prioritaires:**
        1. Volume et symétrie
        2. Insertion des trayons
        3. Veines mammaires
        4. Tempérament calme
        """)
    
    # Statistiques rapides
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE statut = 'active'")
    actives = cursor.fetchone()[0]
    
    st.markdown("### 📊 EN DIRECT")
    st.metric("🐑 Actives", actives)

# ========== PIED DE PAGE ==========
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - RACES ALGÉRIENNES</strong> | Version 4.0</p>
    <p>© 2024 - Tous droits réservés</p>
</div>
""", unsafe_allow_html=True)
