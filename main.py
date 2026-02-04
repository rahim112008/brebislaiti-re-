"""
OVIN MANAGER PRO - Version Complète avec Scanner 3D
Base de données simulée de 50 brebis Hamra Algérienne
"""

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
import cv2

# ========== CONFIGURATION ==========

st.set_page_config(
    page_title="Ovin Manager Pro - Scanner 3D",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personnalisé
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
    .hamra-card {
        background: linear-gradient(135deg, #8B0000 0%, #FF4500 100%);
        color: white;
        padding: 20px;
        border-radius: 15px;
        margin: 15px 0;
        box-shadow: 0 10px 20px rgba(139,0,0,0.2);
    }
    .scanner-card {
        background: linear-gradient(135deg, #1a237e 0%, #283593 100%);
        color: white;
        padding: 20px;
        border-radius: 15px;
        margin: 15px 0;
        box-shadow: 0 10px 20px rgba(26,35,126,0.2);
    }
    .metric-hamra {
        background: linear-gradient(135deg, #FFF5F5 0%, #FFE4E1 100%);
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 5px 15px rgba(139,0,0,0.1);
        border-left: 5px solid #8B0000;
    }
    .scanner-view {
        background: black;
        border-radius: 10px;
        padding: 10px;
        margin: 10px 0;
        text-align: center;
    }
    .scan-progress {
        height: 20px;
        background: linear-gradient(90deg, #1a237e, #283593);
        border-radius: 10px;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

# ========== CRÉATION BASE DE DONNÉES SIMULÉE ==========

def creer_base_hamra():
    """Crée une base de données avec 50 brebis Hamra simulées"""
    conn = sqlite3.connect('hamra_algerienne.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table brebis
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE NOT NULL,
            nom TEXT,
            race TEXT DEFAULT 'HAMRA',
            sous_race TEXT,
            sexe TEXT CHECK(sexe IN ('F', 'M')),
            date_naissance DATE,
            age_mois INTEGER,
            poids FLOAT,
            score_condition INTEGER CHECK(score_condition BETWEEN 1 AND 5),
            couleur_robe TEXT,
            intensite_roux INTEGER CHECK(intensite_roux BETWEEN 1 AND 10),
            cornes BOOLEAN,
            taille_cornes_cm FLOAT,
            longueur_corps_cm FLOAT,
            hauteur_garrot_cm FLOAT,
            largeur_bassin_cm FLOAT,
            tour_poitrine_cm FLOAT,
            notes TEXT,
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
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Table gestations
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_saillie DATE,
            date_mise_bas_prevu DATE,
            nombre_agneaux_prevus INTEGER,
            statut TEXT DEFAULT 'en_cours',
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
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
            fichier_3d_path TEXT,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
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
            date_analyse DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis(id)
        )
    ''')
    
    # Vérifier si la base est vide et la peupler
    cursor.execute("SELECT COUNT(*) FROM brebis")
    count = cursor.fetchone()[0]
    
    if count == 0:
        peupler_base_hamra(cursor, conn)
    
    conn.commit()
    return conn

def peupler_base_hamra(cursor, conn):
    """Peuple la base avec 50 brebis Hamra simulées"""
    
    # Noms traditionnels algériens pour brebis
    noms_femelles = [
        "Lalla", "Zina", "Nour", "Yasmina", "Fatima", "Khadija", "Aicha", "Samira", 
        "Leila", "Soraya", "Nadia", "Rym", "Salima", "Djamila", "Zahra", "Malika",
        "Halima", "Farida", "Rachida", "Safia", "Hind", "Mouna", "Widad", "Nawal",
        "Saida", "Zahia", "Houria", "Taklit", "Aziza", "Rania"
    ]
    
    noms_males = [
        "Sultan", "Amir", "Karim", "Rachid", "Yacine", "Noureddine", "Mohamed", 
        "Ali", "Omar", "Hassan", "Hussein", "Bilal", "Mokhtar", "Salim", "Aziz",
        "Farid", "Kader", "Said", "Tahar", "Mustapha"
    ]
    
    # Créer 35 femelles et 15 mâles
    brebis_data = []
    
    # Femelles
    for i in range(1, 36):
        identifiant = f"HAM-F-2023-{i:03d}"
        nom = random.choice(noms_femelles) + f" {i}"
        age_mois = random.randint(12, 84)  # 1 à 7 ans
        date_naissance = date.today() - timedelta(days=age_mois*30)
        poids = random.uniform(45.0, 65.0)
        score_condition = random.randint(2, 4)
        intensite_roux = random.randint(6, 10)  # Hamra = rousse
        cornes = random.choice([True, False])
        taille_cornes = random.uniform(0, 25) if cornes else 0
        
        # Mensurations
        longueur_corps = random.uniform(90.0, 120.0)
        hauteur_garrot = random.uniform(65.0, 80.0)
        largeur_bassin = random.uniform(35.0, 50.0)
        tour_poitrine = random.uniform(95.0, 120.0)
        
        # Mère et père (simulés)
        mere_id = f"HAM-F-2020-{random.randint(1, 20):03d}" if i > 10 else None
        pere_id = f"HAM-M-2019-{random.randint(1, 10):03d}" if i > 10 else None
        
        brebis_data.append((
            identifiant, nom, 'HAMRA', 'Type El Oued', 'F',
            date_naissance.isoformat(), age_mois, poids, score_condition,
            'Rousse', intensite_roux, cornes, taille_cornes,
            longueur_corps, hauteur_garrot, largeur_bassin, tour_poitrine,
            f"Brebis Hamra {i} - Élevage traditionnel algérien",
            mere_id, pere_id, random.uniform(0.0, 0.15), 'active'
        ))
    
    # Mâles
    for i in range(1, 16):
        identifiant = f"HAM-M-2023-{i:03d}"
        nom = random.choice(noms_males) + f" {i}"
        age_mois = random.randint(12, 72)
        date_naissance = date.today() - timedelta(days=age_mois*30)
        poids = random.uniform(65.0, 90.0)  # Mâles plus lourds
        score_condition = random.randint(3, 5)
        intensite_roux = random.randint(7, 10)
        cornes = True  # Les mâles ont généralement des cornes
        taille_cornes = random.uniform(30.0, 60.0)
        
        # Mensurations (mâles plus grands)
        longueur_corps = random.uniform(100.0, 130.0)
        hauteur_garrot = random.uniform(75.0, 90.0)
        largeur_bassin = random.uniform(40.0, 55.0)
        tour_poitrine = random.uniform(110.0, 135.0)
        
        brebis_data.append((
            identifiant, nom, 'HAMRA', 'Type Ouargla', 'M',
            date_naissance.isoformat(), age_mois, poids, score_condition,
            'Rousse foncée', intensite_roux, cornes, taille_cornes,
            longueur_corps, hauteur_garrot, largeur_bassin, tour_poitrine,
            f"Bélier Hamra {i} - Reproducteur",
            None, None, random.uniform(0.0, 0.05), 'active'
        ))
    
    # Insérer les brebis
    cursor.executemany('''
        INSERT INTO brebis (
            identifiant, nom, race, sous_race, sexe, date_naissance, age_mois, 
            poids, score_condition, couleur_robe, intensite_roux, cornes, 
            taille_cornes_cm, longueur_corps_cm, hauteur_garrot_cm, 
            largeur_bassin_cm, tour_poitrine_cm, notes, mere_id, pere_id, 
            coefficient_consanguinite, statut
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', brebis_data)
    
    # Peupler la production laitière pour les femelles
    production_data = []
    for i in range(1, 36):  # Pour chaque femelle
        brebis_id = i
        # Générer des mesures sur les 90 derniers jours
        for j in range(random.randint(5, 15)):
            date_mesure = date.today() - timedelta(days=random.randint(1, 90))
            quantite = random.uniform(1.5, 3.5)  # Production typique
            taux_mg = random.uniform(6.0, 8.5)   # Hamra: lait riche
            taux_proteine = random.uniform(5.0, 6.5)
            cellules = random.randint(100000, 400000)
            
            production_data.append((
                brebis_id, date_mesure.isoformat(), quantite, taux_mg, 
                taux_proteine, cellules, f"Mesure {j+1}"
            ))
    
    cursor.executemany('''
        INSERT INTO production_lait (
            brebis_id, date_mesure, quantite_litre, taux_matiere_grasse,
            taux_proteine, cellules_somatiques, notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', production_data)
    
    # Peupler les gestations
    gestation_data = []
    femelles_gestantes = random.sample(range(1, 36), 12)  # 12 femelles gestantes
    
    for idx, brebis_id in enumerate(femelles_gestantes):
        date_saillie = date.today() - timedelta(days=random.randint(30, 120))
        date_mise_bas = date_saillie + timedelta(days=150)
        nombre_agneaux = random.choices([1, 2, 3], weights=[0.4, 0.5, 0.1])[0]
        
        gestation_data.append((
            brebis_id, date_saillie.isoformat(), date_mise_bas.isoformat(),
            nombre_agneaux, 'en_cours', f"Gestation {idx+1}"
        ))
    
    cursor.executemany('''
        INSERT INTO gestations (
            brebis_id, date_saillie, date_mise_bas_prevu, nombre_agneaux_prevus,
            statut, notes
        ) VALUES (?, ?, ?, ?, ?, ?)
    ''', gestation_data)
    
    # Peupler les scans 3D
    scans_data = []
    for brebis_id in random.sample(range(1, 51), 20):  # 20 brebis scannées
        date_scan = date.today() - timedelta(days=random.randint(1, 60))
        
        # Générer des points 3D simulés
        points_3d = []
        for _ in range(100):
            points_3d.append({
                'x': random.uniform(-50, 50),
                'y': random.uniform(-30, 30),
                'z': random.uniform(-20, 20),
                'intensity': random.randint(50, 255)
            })
        
        # Mesures estimées à partir du scan
        mesures = {
            'volume': random.uniform(0.08, 0.15),  # en m³
            'surface': random.uniform(1.5, 2.5),    # en m²
            'indice_corporel': random.uniform(2.8, 3.5),
            'symetrie': random.uniform(85, 98)
        }
        
        scans_data.append((
            brebis_id, date_scan.isoformat(), 'Laser 3D',
            json.dumps(points_3d), json.dumps(mesures),
            mesures['volume'], mesures['surface'],
            random.randint(70, 95), f'/scans/brebis_{brebis_id}.ply',
            f'Scan complet - Brebis {brebis_id}'
        ))
    
    cursor.executemany('''
        INSERT INTO scans_3d (
            brebis_id, date_scan, mode_scan, points_3d_json, mesures_json,
            volume_estime, surface_estimee, qualite_scan, fichier_3d_path, notes
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', scans_data)
    
    conn.commit()
    print(f"✅ Base de données peuplée avec {len(brebis_data)} brebis Hamra")

# Initialiser la base
conn = creer_base_hamra()

# ========== MODULE SCANNER 3D ==========

class Scanner3D:
    """Simulateur de scanner 3D pour ovins"""
    
    @staticmethod
    def generer_photo_simulee(brebis_info):
        """Génère une photo simulée d'une brebis"""
        # Créer une image de base
        width, height = 400, 300
        image = Image.new('RGB', (width, height), color='white')
        draw = ImageDraw.Draw(image)
        
        # Couleur selon la race
        if brebis_info.get('race') == 'HAMRA':
            corps_color = (139, 0, 0)  # Rouge foncé pour Hamra
        else:
            corps_color = (255, 255, 255)  # Blanc pour autres races
        
        # Dessiner le corps
        draw.ellipse([100, 80, 300, 200], fill=corps_color, outline='black', width=2)
        
        # Tête
        draw.ellipse([280, 100, 350, 160], fill=corps_color, outline='black', width=2)
        
        # Pattes
        for x in [130, 170, 230, 270]:
            draw.rectangle([x, 200, x+20, 280], fill='black')
        
        # Cornes si présentes
        if brebis_info.get('cornes', False):
            draw.arc([300, 70, 350, 120], start=0, end=180, fill='gray', width=5)
            draw.arc([320, 70, 370, 120], start=0, end=180, fill='gray', width=5)
        
        # Ajouter l'identifiant
        draw.text((10, 10), f"ID: {brebis_info.get('identifiant', 'N/A')}", fill='black')
        draw.text((10, 30), f"Race: {brebis_info.get('race', 'N/A')}", fill='black')
        
        return image
    
    @staticmethod
    def simuler_scan_3d(brebis_info):
        """Simule un scan 3D et retourne des points 3D"""
        np.random.seed(hash(brebis_info.get('identifiant', '')) % 10000)
        
        n_points = 500
        points = []
        
        # Corps ellipsoïde
        for _ in range(n_points):
            # Coordonnées sphériques
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.random.uniform(0, np.pi)
            
            # Rayons approximatifs d'une brebis (en cm)
            rx = brebis_info.get('largeur_bassin_cm', 40) / 2
            ry = brebis_info.get('longueur_corps_cm', 110) / 2
            rz = brebis_info.get('hauteur_garrot_cm', 75) / 2
            
            # Conversion en coordonnées cartésiennes
            x = rx * np.sin(phi) * np.cos(theta) + np.random.normal(0, 2)
            y = ry * np.sin(phi) * np.sin(theta) + np.random.normal(0, 3)
            z = rz * np.cos(phi) + np.random.normal(0, 1.5)
            
            # Ajouter du bruit pour réalisme
            x += np.random.normal(0, 1)
            y += np.random.normal(0, 1.5)
            z += np.random.normal(0, 1)
            
            points.append({
                'x': float(x),
                'y': float(y),
                'z': float(z),
                'intensity': int(np.random.uniform(50, 255))
            })
        
        return points
    
    @staticmethod
    def analyser_points_3d(points):
        """Analyse les points 3D pour extraire des mesures"""
        if not points:
            return {}
        
        # Convertir en arrays numpy
        x = np.array([p['x'] for p in points])
        y = np.array([p['y'] for p in points])
        z = np.array([p['z'] for p in points])
        
        # Calculer les mesures
        longueur = np.max(y) - np.min(y)
        largeur = np.max(x) - np.min(x)
        hauteur = np.max(z) - np.min(z)
        
        # Volume approximatif (ellipsoïde)
        volume = (4/3) * np.pi * (longueur/2) * (largeur/2) * (hauteur/2) / 1000000  # en m³
        
        # Surface approximative
        a, b, c = longueur/2, largeur/2, hauteur/2
        surface = 4 * np.pi * (a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3)**(1/1.6) / 10000  # en m²
        
        # Indice corporel
        indice_corporel = (longueur * largeur * hauteur) ** (1/3)
        
        # Symétrie (comparaison gauche/droite)
        symetrie = 100 - np.abs(np.mean(x[x > 0]) + np.mean(x[x < 0])) * 10
        
        return {
            'longueur_cm': round(float(longueur), 1),
            'largeur_cm': round(float(largeur), 1),
            'hauteur_cm': round(float(hauteur), 1),
            'volume_m3': round(float(volume), 4),
            'surface_m2': round(float(surface), 2),
            'indice_corporel': round(float(indice_corporel), 2),
            'symetrie_percent': round(float(symetrie), 1)
        }
    
    @staticmethod
    def creer_visualisation_3d(points):
        """Crée une visualisation 3D des points"""
        if len(points) < 10:
            return None
        
        # Extraire les coordonnées
        x = [p['x'] for p in points]
        y = [p['y'] for p in points]
        z = [p['z'] for p in points]
        intensities = [p['intensity'] for p in points]
        
        # Créer le graphique 3D
        fig = go.Figure(data=[go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=3,
                color=intensities,
                colorscale='Viridis',
                opacity=0.8,
                showscale=True
            )
        )])
        
        fig.update_layout(
            title="Scan 3D - Reconstruction",
            scene=dict(
                xaxis_title='Largeur (cm)',
                yaxis_title='Longueur (cm)',
                zaxis_title='Hauteur (cm)',
                aspectmode='data'
            ),
            width=800,
            height=600
        )
        
        return fig
    
    @staticmethod
    def generer_mesh_3d(points):
        """Génère un maillage 3D simplifié à partir des points"""
        # Pour la démo, retourner une sphère
        u = np.linspace(0, 2*np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        
        x = 20 * np.outer(np.cos(u), np.sin(v))
        y = 40 * np.outer(np.sin(u), np.sin(v))
        z = 15 * np.outer(np.ones(np.size(u)), np.cos(v))
        
        return x, y, z

# ========== PAGES DE L'APPLICATION ==========

def page_accueil():
    """Page d'accueil avec vue d'ensemble"""
    st.markdown('<h1 class="main-header">🐑 OVIN MANAGER PRO - HAMRA ALGÉRIENNE</h1>', unsafe_allow_html=True)
    st.markdown("**Système de gestion et d'analyse scientifique du troupeau Hamra Algérien**")
    
    # Métriques principales
    cursor = conn.cursor()
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-hamra'>
            <h3>🐑 TOTAL BREBIS</h3>
            <h2>{total}</h2>
            <p>Race Hamra Algérienne</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F'")
        femelles = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-hamra'>
            <h3>♀️ FEMELLES</h3>
            <h2>{femelles}</h2>
            <p>Productrices</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        gestations = cursor.fetchone()[0]
        st.markdown(f"""
        <div class='metric-hamra'>
            <h3>🤰 GESTATIONS</h3>
            <h2>{gestations}</h2>
            <p>En cours</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        cursor.execute("SELECT AVG(quantite_litre) FROM production_lait WHERE date_mesure > date('now', '-30 days')")
        prod = cursor.fetchone()[0] or 0
        st.markdown(f"""
        <div class='metric-hamra'>
            <h3>🥛 PRODUCTION</h3>
            <h2>{prod:.1f} L/j</h2>
            <p>Moyenne 30 jours</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Graphiques
    st.markdown("### 📊 ANALYSE DU TROUPEAU HAMRA")
    
    col_chart1, col_chart2 = st.columns(2)
    
    with col_chart1:
        # Répartition par âge
        cursor.execute("""
            SELECT 
                CASE 
                    WHEN age_mois < 12 THEN '0-1 an'
                    WHEN age_mois < 24 THEN '1-2 ans'
                    WHEN age_mois < 36 THEN '2-3 ans'
                    WHEN age_mois < 60 THEN '3-5 ans'
                    ELSE '5+ ans'
                END as tranche_age,
                COUNT(*) as count
            FROM brebis
            GROUP BY tranche_age
            ORDER BY tranche_age
        """)
        data_age = cursor.fetchall()
        
        if data_age:
            df_age = pd.DataFrame(data_age, columns=['Tranche', 'Nombre'])
            fig = px.pie(df_age, values='Nombre', names='Tranche', 
                        title="Pyramide des âges", hole=0.4)
            st.plotly_chart(fig, use_container_width=True)
    
    with col_chart2:
        # Production laitière moyenne
        cursor.execute("""
            SELECT strftime('%Y-%m', date_mesure) as mois,
                   AVG(quantite_litre) as moyenne_lait,
                   AVG(taux_matiere_grasse) as moyenne_mg
            FROM production_lait
            WHERE date_mesure > date('now', '-6 months')
            GROUP BY mois
            ORDER BY mois
        """)
        data_prod = cursor.fetchall()
        
        if data_prod:
            df_prod = pd.DataFrame(data_prod, columns=['Mois', 'Lait (L)', 'MG (%)'])
            fig = px.line(df_prod, x='Mois', y=['Lait (L)', 'MG (%)'],
                         title="Évolution production 6 mois", markers=True)
            st.plotly_chart(fig, use_container_width=True)
    
    # Dernières alertes
    st.markdown("### ⚠️ ALERTES ET SUIVI")
    
    # Gestations proches
    cursor.execute("""
        SELECT b.identifiant, b.nom, g.date_mise_bas_prevu,
               julianday(g.date_mise_bas_prevu) - julianday('now') as jours_restants
        FROM gestations g
        JOIN brebis b ON g.brebis_id = b.id
        WHERE g.statut = 'en_cours'
        AND jours_restants BETWEEN 0 AND 14
        ORDER BY jours_restants
    """)
    
    alertes = cursor.fetchall()
    
    if alertes:
        for alerte in alertes:
            jours = int(alerte[3])
            if jours <= 7:
                st.warning(f"🚨 **{alerte[1]}** ({alerte[0]}) - Mise bas dans {jours} jours!")
            else:
                st.info(f"📅 **{alerte[1]}** ({alerte[0]}) - Mise bas dans {jours} jours")
    else:
        st.success("✅ Aucune mise bas imminente")

def page_scanner_3d():
    """Page du scanner 3D"""
    st.markdown('<h2 class="section-header">📐 SCANNER 3D - MORPHOMÉTRIE AVANCÉE</h2>', unsafe_allow_html=True)
    
    # Sélection de la brebis
    cursor = conn.cursor()
    cursor.execute("SELECT id, identifiant, nom, race FROM brebis ORDER BY nom")
    brebis_list = cursor.fetchall()
    
    if not brebis_list:
        st.warning("Aucune brebis dans la base de données")
        return
    
    col_sel1, col_sel2 = st.columns([2, 1])
    
    with col_sel1:
        brebis_option = st.selectbox(
            "SÉLECTIONNEZ UNE BREBIS À SCANNER:",
            [f"{b[2]} ({b[1]}) - {b[3]}" for b in brebis_list]
        )
    
    with col_sel2:
        mode_scan = st.selectbox(
            "MODE DE SCAN:",
            ["Laser haute précision", "Photogrammétrie", "Scanner portable", "Simulation"]
        )
    
    if brebis_option:
        # Extraire l'ID
        brebis_id = int(brebis_option.split('(')[1].split(')')[0].split('-')[-1])
        
        # Récupérer les infos de la brebis
        cursor.execute("SELECT * FROM brebis WHERE id = ?", (brebis_id,))
        brebis_info = cursor.fetchone()
        columns = [desc[0] for desc in cursor.description]
        brebis_dict = dict(zip(columns, brebis_info))
        
        # Onglets du scanner
        scan_tabs = st.tabs(["📸 PHOTO", "🎯 SCAN 3D", "📏 MESURES", "📊 ANALYSE"])
        
        with scan_tabs[0]:
            st.markdown("### 📸 PHOTOGRAPHIE DE LA BREBIS")
            
            # Générer une photo simulée
            photo = Scanner3D.generer_photo_simulee(brebis_dict)
            
            col_photo1, col_photo2 = st.columns([2, 1])
            
            with col_photo1:
                st.image(photo, caption=f"Photo simulée - {brebis_dict['nom']}", use_column_width=True)
                
                # Boutons photo
                col_btn1, col_btn2 = st.columns(2)
                with col_btn1:
                    if st.button("📸 Prendre une photo", type="primary"):
                        st.success("Photo prise et sauvegardée!")
                
                with col_btn2:
                    if st.button("🔄 Régénérer"):
                        st.rerun()
            
            with col_photo2:
                st.markdown("""
                <div class='scanner-card'>
                    <h4>📷 INFORMATIONS PHOTO</h4>
                    <p><strong>Résolution:</strong> 12 MP</p>
                    <p><strong>ISO:</strong> 400</p>
                    <p><strong>Ouverture:</strong> f/5.6</p>
                    <p><strong>Distance:</strong> 3 m</p>
                    <p><strong>Éclairage:</strong> Naturel</p>
                </div>
                """, unsafe_allow_html=True)
                
                # Télécharger la photo
                buffered = io.BytesIO()
                photo.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                
                st.markdown(f"""
                <a href="data:image/png;base64,{img_str}" download="brebis_{brebis_dict['identifiant']}.png">
                    <button style='background: #8B0000; color: white; padding: 10px 20px; border: none; border-radius: 5px; cursor: pointer; width: 100%;'>
                        📥 Télécharger la photo
                    </button>
                </a>
                """, unsafe_allow_html=True)
        
        with scan_tabs[1]:
            st.markdown("### 🎯 SCAN 3D EN TEMPS RÉEL")
            
            # Simulation de progression du scan
            scan_progress = st.slider("Progression du scan:", 0, 100, 50)
            
            col_scan1, col_scan2 = st.columns([3, 1])
            
            with col_scan1:
                # Zone de visualisation du scan
                st.markdown('<div class="scanner-view">', unsafe_allow_html=True)
                
                # Barre de progression stylisée
                st.markdown(f"""
                <div style="width: 100%; background: #f0f0f0; border-radius: 10px; margin: 20px 0;">
                    <div class="scan-progress" style="width: {scan_progress}%;"></div>
                </div>
                <p style="color: white; font-size: 1.2em;">Scan en cours: {scan_progress}%</p>
                """, unsafe_allow_html=True)
                
                # Points de scan simulés
                if scan_progress > 0:
                    points = Scanner3D.simuler_scan_3d(brebis_dict)
                    
                    # Afficher un aperçu des premiers points
                    st.markdown("**Points 3D capturés:**")
                    df_points = pd.DataFrame(points[:10])
                    st.dataframe(df_points)
                    
                    # Bouton de scan
                    if st.button("🚀 Démarrer le scan 3D", type="primary"):
                        with st.spinner("Scan en cours... Cette opération peut prendre 30 secondes"):
                            # Simulation du scan
                            for i in range(1, 101):
                                st.progress(i/100)
                                # Petite pause pour l'effet visuel
                                if i % 10 == 0:
                                    st.write(f"Capture {i}%...")
                            
                            # Sauvegarder le scan
                            cursor.execute('''
                                INSERT INTO scans_3d (brebis_id, date_scan, mode_scan, points_3d_json, qualite_scan, notes)
                                VALUES (?, ?, ?, ?, ?, ?)
                            ''', (
                                brebis_id,
                                date.today().isoformat(),
                                mode_scan,
                                json.dumps(points[:100]),  # Limité à 100 points pour la démo
                                85,
                                f"Scan {mode_scan} - {brebis_dict['nom']}"
                            ))
                            conn.commit()
                            st.success("✅ Scan 3D terminé et sauvegardé!")
                
                st.markdown('</div>', unsafe_allow_html=True)
            
            with col_scan2:
                st.markdown("""
                <div class='scanner-card'>
                    <h4>⚙️ PARAMÈTRES SCAN</h4>
                    <p><strong>Mode:</strong> {}</p>
                    <p><strong>Précision:</strong> ±0.5 mm</p>
                    <p><strong>Points/s:</strong> 50,000</p>
                    <p><strong>Portée:</strong> 5 m</p>
                    <p><strong>Laser:</strong> Classe 1</p>
                </div>
                """.format(mode_scan), unsafe_allow_html=True)
        
        with scan_tabs[2]:
            st.markdown("### 📏 MESURES MORPHOMÉTRIQUES")
            
            # Récupérer ou générer des mesures
            cursor.execute("SELECT * FROM scans_3d WHERE brebis_id = ? ORDER BY date_scan DESC LIMIT 1", (brebis_id,))
            scan_data = cursor.fetchone()
            
            if scan_data:
                # Mesures existantes
                scan_columns = [desc[0] for desc in cursor.description]
                scan_dict = dict(zip(scan_columns, scan_data))
                
                mesures = json.loads(scan_dict['mesures_json']) if scan_dict['mesures_json'] else {}
            else:
                # Générer de nouvelles mesures
                points = Scanner3D.simuler_scan_3d(brebis_dict)
                mesures = Scanner3D.analyser_points_3d(points)
            
            # Afficher les mesures
            if mesures:
                col_mes1, col_mes2, col_mes3 = st.columns(3)
                
                with col_mes1:
                    st.metric("Longueur", f"{mesures.get('longueur_cm', 0):.1f} cm")
                    st.metric("Volume", f"{mesures.get('volume_m3', 0):.3f} m³")
                
                with col_mes2:
                    st.metric("Largeur", f"{mesures.get('largeur_cm', 0):.1f} cm")
                    st.metric("Surface", f"{mesures.get('surface_m2', 0):.2f} m²")
                
                with col_mes3:
                    st.metric("Hauteur", f"{mesures.get('hauteur_cm', 0):.1f} cm")
                    st.metric("Symétrie", f"{mesures.get('symetrie_percent', 0):.1f}%")
                
                # Graphique des proportions
                fig = go.Figure(data=[go.Bar(
                    x=['Longueur', 'Largeur', 'Hauteur'],
                    y=[mesures.get('longueur_cm', 0), 
                       mesures.get('largeur_cm', 0), 
                       mesures.get('hauteur_cm', 0)],
                    marker_color=['#8B0000', '#FF4500', '#FF8C00']
                )])
                
                fig.update_layout(
                    title="Proportions corporelles",
                    yaxis_title="Centimètres (cm)"
                )
                
                st.plotly_chart(fig, use_container_width=True)
            
            # Comparaison avec les standards Hamra
            st.markdown("### 🎯 COMPARAISON AVEC STANDARDS HAMRA")
            
            standards_hamra = {
                'longueur_cm': {'min': 95, 'max': 125, 'optimum': 110},
                'largeur_cm': {'min': 35, 'max': 55, 'optimum': 45},
                'hauteur_cm': {'min': 65, 'max': 85, 'optimum': 75},
                'indice_corporel': {'min': 2.5, 'max': 3.5, 'optimum': 3.0}
            }
            
            df_compare = pd.DataFrame([
                {
                    'Mesure': 'Longueur',
                    'Valeur': mesures.get('longueur_cm', 0),
                    'Standard Min': standards_hamra['longueur_cm']['min'],
                    'Standard Max': standards_hamra['longueur_cm']['max'],
                    'Optimum': standards_hamra['longueur_cm']['optimum']
                },
                {
                    'Mesure': 'Largeur',
                    'Valeur': mesures.get('largeur_cm', 0),
                    'Standard Min': standards_hamra['largeur_cm']['min'],
                    'Standard Max': standards_hamra['largeur_cm']['max'],
                    'Optimum': standards_hamra['largeur_cm']['optimum']
                },
                {
                    'Mesure': 'Hauteur',
                    'Valeur': mesures.get('hauteur_cm', 0),
                    'Standard Min': standards_hamra['hauteur_cm']['min'],
                    'Standard Max': standards_hamra['hauteur_cm']['max'],
                    'Optimum': standards_hamra['hauteur_cm']['optimum']
                }
            ])
            
            # Graphique de comparaison
            fig = go.Figure()
            
            for idx, row in df_compare.iterrows():
                fig.add_trace(go.Scatter(
                    x=[row['Mesure']],
                    y=[row['Valeur']],
                    mode='markers',
                    name='Mesure',
                    marker=dict(size=15, color='#8B0000')
                ))
                
                # Barre d'intervalle optimal
                fig.add_trace(go.Scatter(
                    x=[row['Mesure'], row['Mesure']],
                    y=[row['Standard Min'], row['Standard Max']],
                    mode='lines',
                    name='Intervalle normal',
                    line=dict(color='green', width=10),
                    showlegend=(idx == 0)
                ))
            
            fig.update_layout(
                title="Comparaison avec standards race Hamra",
                yaxis_title="Centimètres (cm)",
                showlegend=True
            )
            
            st.plotly_chart(fig, use_container_width=True)
        
        with scan_tabs[3]:
            st.markdown("### 📊 ANALYSE AVANCÉE 3D")
            
            # Générer une visualisation 3D
            points = Scanner3D.simuler_scan_3d(brebis_dict)
            fig_3d = Scanner3D.creer_visualisation_3d(points)
            
            if fig_3d:
                st.plotly_chart(fig_3d, use_container_width=True)
            
            # Statistiques du scan
            st.markdown("#### 📈 STATISTIQUES DU NUAGE DE POINTS")
            
            if points:
                stats = {
                    'Nombre de points': len(points),
                    'Densité (pts/cm³)': round(len(points) / 1000, 1),
                    'Précision estimée': '±0.8 mm',
                    'Couverture surface': '92%',
                    'Points manquants': '8%'
                }
                
                col_stats1, col_stats2 = st.columns(2)
                
                with col_stats1:
                    for key, value in list(stats.items())[:3]:
                        st.metric(key, value)
                
                with col_stats2:
                    for key, value in list(stats.items())[3:]:
                        st.metric(key, value)
                
                # Analyse de forme
                st.markdown("#### 🎯 ANALYSE DE FORME ET CONFORMATION")
                
                # Générer un maillage 3D simplifié
                x, y, z = Scanner3D.generer_mesh_3d(points)
                
                fig_mesh = go.Figure(data=[go.Surface(
                    x=x, y=y, z=z,
                    colorscale='Reds',
                    opacity=0.8,
                    contours={
                        "z": {"show": True, "usecolormap": True, "highlightcolor": "limegreen", "project": {"z": True}}
                    }
                )])
                
                fig_mesh.update_layout(
                    title="Modélisation 3D - Surface reconstruite",
                    scene=dict(
                        xaxis_title='Largeur',
                        yaxis_title='Longueur',
                        zaxis_title='Hauteur'
                    ),
                    width=800,
                    height=500
                )
                
                st.plotly_chart(fig_mesh, use_container_width=True)
                
                # Recommandations
                st.markdown("#### 💡 RECOMMANDATIONS")
                
                recommendations = []
                if mesures.get('indice_corporel', 0) < 2.8:
                    recommendations.append("Améliorer l'alimentation pour augmenter l'indice corporel")
                if mesures.get('symetrie_percent', 100) < 90:
                    recommendations.append("Surveiller la symétrie - possible problème locomoteur")
                if len(points) < 300:
                    recommendations.append("Qualité de scan insuffisante - recommandé de refaire le scan")
                
                if recommendations:
                    for rec in recommendations:
                        st.warning(f"⚠️ {rec}")
                else:
                    st.success("✅ Conformation optimale - Animal bien conformé selon les standards Hamra")

def page_gestion():
    """Page de gestion du troupeau"""
    st.markdown('<h2 class="section-header">📊 GESTION DU TROUPEAU HAMRA</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🐑 LISTE COMPLÈTE", "📈 STATISTIQUES", "🔍 RECHERCHE", "📤 EXPORT"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT identifiant, nom, sexe, age_mois, poids, score_condition, 
                   intensite_roux, cornes, statut, date_naissance
            FROM brebis
            ORDER BY identifiant
        """)
        
        brebis_data = cursor.fetchall()
        columns = ['ID', 'Nom', 'Sexe', 'Âge (mois)', 'Poids (kg)', 'Score', 
                  'Intensité Roux', 'Cornes', 'Statut', 'Naissance']
        
        df = pd.DataFrame(brebis_data, columns=columns)
        
        # Recherche
        recherche = st.text_input("🔍 Rechercher une brebis:")
        if recherche:
            df = df[df.apply(lambda row: row.astype(str).str.contains(recherche, case=False).any(), axis=1)]
        
        st.dataframe(df, use_container_width=True, height=600)
        
        st.metric("Brebis affichées", len(df))
    
    with tab2:
        st.markdown("### 📊 ANALYSE STATISTIQUE DU TROUPEAU")
        
        col_stat1, col_stat2, col_stat3 = st.columns(3)
        
        with col_stat1:
            cursor.execute("SELECT AVG(poids) FROM brebis WHERE sexe = 'F'")
            poids_f = cursor.fetchone()[0]
            st.metric("Poids moyen ♀️", f"{poids_f:.1f} kg")
        
        with col_stat2:
            cursor.execute("SELECT AVG(poids) FROM brebis WHERE sexe = 'M'")
            poids_m = cursor.fetchone()[0]
            st.metric("Poids moyen ♂️", f"{poids_m:.1f} kg")
        
        with col_stat3:
            cursor.execute("SELECT AVG(intensite_roux) FROM brebis")
            roux = cursor.fetchone()[0]
            st.metric("Intensité roux", f"{roux:.1f}/10")
        
        # Graphique de distribution des poids
        cursor.execute("SELECT poids, sexe FROM brebis")
        poids_data = cursor.fetchall()
        
        if poids_data:
            df_poids = pd.DataFrame(poids_data, columns=['Poids', 'Sexe'])
            
            fig = px.histogram(df_poids, x='Poids', color='Sexe',
                             title="Distribution des poids par sexe",
                             nbins=20,
                             color_discrete_map={'F': '#FF6B6B', 'M': '#4ECDC4'})
            
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.markdown("### 🔍 RECHERCHE AVANCÉE")
        
        col_search1, col_search2 = st.columns(2)
        
        with col_search1:
            min_poids = st.number_input("Poids minimum (kg)", 0, 200, 40)
            max_poids = st.number_input("Poids maximum (kg)", 0, 200, 100)
            
            score_min = st.slider("Score condition minimum", 1, 5, 2)
        
        with col_search2:
            intensite_min = st.slider("Intensité rousse minimum", 1, 10, 5)
            avec_cornes = st.selectbox("Cornes", ["Tous", "Avec cornes", "Sans cornes"])
            statut = st.selectbox("Statut", ["Tous", "active", "retrait", "malade"])
        
        # Construire la requête
        query = "SELECT * FROM brebis WHERE 1=1"
        params = []
        
        query += " AND poids BETWEEN ? AND ?"
        params.extend([min_poids, max_poids])
        
        query += " AND score_condition >= ?"
        params.append(score_min)
        
        query += " AND intensite_roux >= ?"
        params.append(intensite_min)
        
        if avec_cornes == "Avec cornes":
            query += " AND cornes = 1"
        elif avec_cornes == "Sans cornes":
            query += " AND cornes = 0"
        
        if statut != "Tous":
            query += " AND statut = ?"
            params.append(statut)
        
        cursor.execute(query, params)
        resultats = cursor.fetchall()
        
        st.metric("Résultats trouvés", len(resultats))
        
        if resultats:
            df_result = pd.DataFrame(resultats, columns=[desc[0] for desc in cursor.description])
            st.dataframe(df_result[['identifiant', 'nom', 'sexe', 'poids', 'score_condition', 'intensite_roux']])
    
    with tab4:
        st.markdown("### 📤 EXPORTATION DES DONNÉES")
        
        export_format = st.selectbox("Format d'export", ["CSV", "Excel", "JSON", "PDF"])
        
        if st.button("📥 Générer l'export", type="primary"):
            cursor.execute("SELECT * FROM brebis")
            data = cursor.fetchall()
            columns = [desc[0] for desc in cursor.description]
            df_export = pd.DataFrame(data, columns=columns)
            
            if export_format == "CSV":
                csv = df_export.to_csv(index=False)
                st.download_button(
                    label="Télécharger CSV",
                    data=csv,
                    file_name="troupeau_hamra.csv",
                    mime="text/csv"
                )
            
            elif export_format == "JSON":
                json_data = df_export.to_json(orient='records', indent=2)
                st.download_button(
                    label="Télécharger JSON",
                    data=json_data,
                    file_name="troupeau_hamra.json",
                    mime="application/json"
                )

def page_production():
    """Page de suivi de production"""
    st.markdown('<h2 class="section-header">🥛 SUIVI DE PRODUCTION LAITIÈRE</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📝 SAISIE", "📈 ANALYSE", "🏆 CLASSEMENT"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom FROM brebis WHERE sexe = 'F'")
        femelles = cursor.fetchall()
        
        if femelles:
            with st.form("form_production"):
                brebis_sel = st.selectbox("Sélectionner une brebis", 
                                        [f"{f[1]} (ID:{f[0]})" for f in femelles])
                
                col_prod1, col_prod2, col_prod3 = st.columns(3)
                
                with col_prod1:
                    quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5, 0.1)
                
                with col_prod2:
                    mg = st.number_input("Matière grasse %", 0.0, 20.0, 7.2, 0.1)
                
                with col_prod3:
                    proteine = st.number_input("Protéine %", 0.0, 20.0, 5.5, 0.1)
                
                cellules = st.number_input("Cellules somatiques (x1000)", 0, 1000, 200)
                date_mesure = st.date_input("Date", value=date.today())
                
                if st.form_submit_button("💾 Enregistrer", type="primary"):
                    brebis_id = int(brebis_sel.split("ID:")[1].strip(")"))
                    
                    cursor.execute('''
                        INSERT INTO production_lait 
                        (brebis_id, date_mesure, quantite_litre, taux_matiere_grasse, 
                         taux_proteine, cellules_somatiques)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (brebis_id, date_mesure.isoformat(), quantite, mg, proteine, cellules*1000))
                    conn.commit()
                    st.success("✅ Production enregistrée!")
        else:
            st.warning("Aucune brebis femelle enregistrée")
    
    with tab2:
        cursor.execute("""
            SELECT strftime('%Y-%m', p.date_mesure) as mois,
                   AVG(p.quantite_litre) as lait_moyen,
                   AVG(p.taux_matiere_grasse) as mg_moyen,
                   AVG(p.taux_proteine) as proteine_moyenne,
                   COUNT(*) as nb_mesures
            FROM production_lait p
            WHERE p.date_mesure > date('now', '-12 months')
            GROUP BY mois
            ORDER BY mois
        """)
        
        data_tendance = cursor.fetchall()
        
        if data_tendance:
            df_tendance = pd.DataFrame(data_tendance, 
                                     columns=['Mois', 'Lait (L)', 'MG (%)', 'Protéine (%)', 'Mesures'])
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df_tendance['Mois'], y=df_tendance['Lait (L)'], 
                                   name='Lait (L)', line=dict(color='blue')))
            fig.add_trace(go.Scatter(x=df_tendance['Mois'], y=df_tendance['MG (%)'], 
                                   name='MG (%)', yaxis='y2', line=dict(color='red')))
            
            fig.update_layout(
                title="Évolution de la production sur 12 mois",
                yaxis=dict(title="Lait (L)"),
                yaxis2=dict(title="%", overlaying='y', side='right'),
                hovermode='x unified'
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Analyse qualité
            st.markdown("### 🧪 ANALYSE QUALITÉ DU LAIT")
            
            col_qual1, col_qual2, col_qual3 = st.columns(3)
            
            with col_qual1:
                mg_moy = df_tendance['MG (%)'].mean()
                st.metric("MG moyenne", f"{mg_moy:.1f}%")
            
            with col_qual2:
                prot_moy = df_tendance['Protéine (%)'].mean()
                st.metric("Protéine moyenne", f"{prot_moy:.1f}%")
            
            with col_qual3:
                ratio = mg_moy / prot_moy if prot_moy > 0 else 0
                st.metric("Ratio MG/Protéine", f"{ratio:.2f}")
    
    with tab3:
        st.markdown("### 🏆 CLASSEMENT DES MEILLEURES PRODUCTRICES")
        
        cursor.execute("""
            SELECT b.nom, b.identifiant,
                   AVG(p.quantite_litre) as moyenne_lait,
                   AVG(p.taux_matiere_grasse) as moyenne_mg,
                   AVG(p.taux_proteine) as moyenne_proteine,
                   COUNT(*) as nb_mesures
            FROM production_lait p
            JOIN brebis b ON p.brebis_id = b.id
            WHERE p.date_mesure > date('now', '-90 days')
            GROUP BY b.id
            HAVING nb_mesures >= 5
            ORDER BY moyenne_lait DESC
            LIMIT 10
        """)
        
        top_prod = cursor.fetchall()
        
        if top_prod:
            df_top = pd.DataFrame(top_prod, 
                                columns=['Nom', 'ID', 'Lait moyen (L)', 'MG (%)', 'Protéine (%)', 'Mesures'])
            
            st.dataframe(df_top.style.highlight_max(subset=['Lait moyen (L)'], color='lightgreen'))
            
            # Graphique top 5
            fig = px.bar(df_top.head(5), x='Nom', y='Lait moyen (L)',
                        color='MG (%)',
                        title="Top 5 productrices - 90 derniers jours",
                        hover_data=['MG (%)', 'Protéine (%)'])
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Données insuffisantes pour le classement")

# ========== NAVIGATION PRINCIPALE ==========

# Sidebar
with st.sidebar:
    st.markdown("""
    <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #8B0000 0%, #FF4500 100%); 
                color: white; border-radius: 10px; margin-bottom: 20px;'>
        <h2>🐑 HAMRA ALGÉRIENNE</h2>
        <p>Élevage scientifique</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("### 📍 NAVIGATION")
    
    page = st.radio(
        "MENU PRINCIPAL",
        ["🏠 ACCUEIL", 
         "📐 SCANNER 3D", 
         "📊 GESTION", 
         "🥛 PRODUCTION", 
         "🤰 GESTATION",
         "🧬 GÉNÉTIQUE",
         "⚙️ PARAMÈTRES"]
    )
    
    st.markdown("---")
    
    # Statistiques en temps réel
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F' AND statut = 'active'")
    femelles_actives = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
    gest_actives = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT brebis_id) FROM scans_3d")
    scans_realises = cursor.fetchone()[0]
    
    st.markdown("### 📊 EN DIRECT")
    st.metric("♀️ Actives", femelles_actives)
    st.metric("🤰 Gestations", gest_actives)
    st.metric("📐 Scans 3D", scans_realises)

# Navigation des pages
if page == "🏠 ACCUEIL":
    page_accueil()
elif page == "📐 SCANNER 3D":
    page_scanner_3d()
elif page == "📊 GESTION":
    page_gestion()
elif page == "🥛 PRODUCTION":
    page_production()
elif page == "🤰 GESTATION":
    st.info("Page gestation - À compléter")
elif page == "🧬 GÉNÉTIQUE":
    st.info("Page génétique - À compléter")
elif page == "⚙️ PARAMÈTRES":
    st.info("Page paramètres - À compléter")

# Pied de page
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🐑 <strong>OVIN MANAGER PRO - HAMRA ALGÉRIENNE</strong> | Base simulée: 50 brebis</p>
    <p>📐 Module Scanner 3D inclus | Version 3.0 | © 2024</p>
</div>
""", unsafe_allow_html=True)
