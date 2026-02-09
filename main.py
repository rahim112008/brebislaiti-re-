"""
EXPERT OVIN DZ PRO - VERSION STREAMLIT CLOUD COMPATIBLE
Système de gestion ovine avec mesures morphométriques.
Version allégée sans packages lourds.
"""

import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import os
import json
import random
from datetime import datetime, date
import io
import base64
import math

# ============================================================================
# 1. CONFIGURATION INITIALE
# ============================================================================

# Configuration de la page
st.set_page_config(
    page_title="Expert Ovin DZ",
    page_icon="🐑",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Titre principal
st.title("🐑 EXPERT OVIN DZ PRO")
st.markdown("**Système intégré de gestion ovine algérienne avec mesures morphométriques**")

# ============================================================================
# 2. BASE DE DONNEES
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "ovin_master.db"):
        self.db_path = db_path
        if not os.path.exists('data'): 
            os.makedirs('data', exist_ok=True)
        self.db_path = os.path.join('data', db_path)
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
    
    def execute(self, query: str, params: tuple = ()):
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except Exception as e:
            st.error(f"Erreur SQL: {e}")
            return None
    
    def fetch(self, query: str, params: tuple = ()):
        try:
            return pd.read_sql_query(query, self.conn, params=params)
        except Exception as e:
            st.error(f"Erreur fetch: {e}")
            return pd.DataFrame()

def init_database(db):
    """Initialise les tables de la base de données"""
    tables = [
        # Table principale des animaux
        """CREATE TABLE IF NOT EXISTS animaux (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE NOT NULL,
            nom TEXT,
            race TEXT,
            sexe TEXT CHECK(sexe IN ('Femelle', 'Mâle')),
            age_annees REAL,
            date_naissance DATE,
            date_enregistrement DATE DEFAULT CURRENT_DATE,
            -- Mesures corporelles
            hauteur_garrot_cm REAL,
            longueur_corps_cm REAL,
            tour_poitrine_cm REAL,
            largeur_bassin_cm REAL,
            poids_kg REAL,
            note_condition INTEGER CHECK(note_condition BETWEEN 1 AND 5),
            -- Mesures mammaires (pour femelles)
            largeur_mamelle_cm REAL,
            hauteur_mamelle_cm REAL,
            distance_tetines_cm REAL,
            diametre_tetine_cm REAL,
            symetrie_score INTEGER CHECK(symetrie_score BETWEEN 1 AND 10),
            volume_mamelle_ml REAL,
            score_morpho_mamelle INTEGER CHECK(score_morpho_mamelle BETWEEN 1 AND 9),
            -- Image et détection
            image_path TEXT,
            race_detectee TEXT,
            confiance_detection REAL,
            date_analyse DATE
        )""",
        
        # Table de production laitière
        """CREATE TABLE IF NOT EXISTS production_laitiere (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            animal_id TEXT NOT NULL,
            date_controle DATE,
            quantite_lait_l REAL,
            taux_butyreux REAL,
            taux_proteique REAL,
            cellules_somatiques INTEGER,
            FOREIGN KEY (animal_id) REFERENCES animaux(identifiant)
        )""",
        
        # Table des gestations
        """CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            animal_id TEXT NOT NULL,
            date_saillie DATE,
            date_mise_bas_prevue DATE,
            statut TEXT CHECK(statut IN ('En cours', 'Terminée', 'Échouée')),
            nombre_agneaux INTEGER,
            FOREIGN KEY (animal_id) REFERENCES animaux(identifiant)
        )""",
        
        # Table des races algériennes
        """CREATE TABLE IF NOT EXISTS races_algeriennes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nom TEXT UNIQUE NOT NULL,
            description TEXT,
            region TEXT,
            hauteur_moyenne_cm REAL,
            poids_moyen_kg REAL,
            aptitude TEXT,
            image_exemple TEXT
        )"""
    ]
    
    for table in tables:
        db.execute(table)
    
    # Insérer les races algériennes si la table est vide
    races = db.fetch("SELECT COUNT(*) as count FROM races_algeriennes")
    if races.empty or races['count'][0] == 0:
        races_data = [
            ("Ouled Djellal", "Race la plus répandue, excellente adaptation aux steppes", "Steppes", 75, 70, "Mixte (lait-viande)", "🟤"),
            ("Rembi", "Race à viande de qualité supérieure, croissance rapide", "Hauts plateaux", 70, 65, "Viande", "⚫"),
            ("Hamra", "Race rousse très rustique adaptée aux zones arides", "Sud", 68, 60, "Rustique", "🔴"),
            ("D'man", "Race prolifique des oasis, bonne laitière", "Oasis", 65, 55, "Lait", "🟡"),
            ("Berbère", "Race ancienne très rustique des montagnes", "Kabylie", 72, 68, "Mixte", "⚪"),
            ("Sidaho", "Race de taille moyenne, bonne conformation", "Est", 73, 67, "Viande", "🟠"),
            ("Touareg", "Race du désert, très résistante", "Sahara", 66, 58, "Rustique", "🔵")
        ]
        
        for race in races_data:
            db.execute(
                """INSERT OR IGNORE INTO races_algeriennes 
                   (nom, description, region, hauteur_moyenne_cm, poids_moyen_kg, aptitude, image_exemple) 
                   VALUES (?, ?, ?, ?, ?, ?, ?)""",
                race
            )

def creer_donnees_demo(db):
    """Crée des données de démonstration"""
    # Ajouter quelques animaux de démo
    animaux_demo = [
        ("DZ-2024-001", "Bella", "Ouled Djellal", "Femelle", 3.5, "2021-06-15",
         74.5, 82.3, 95.2, 32.1, 68.5, 4,
         25.4, 18.7, 10.5, 2.1, 8, 1250, 7,
         None, "Ouled Djellal", 0.92, date.today()),
        
        ("DZ-2024-002", "Rocky", "Rembi", "Mâle", 4.0, "2020-03-22",
         71.2, 85.6, 98.7, 35.4, 72.3, 5,
         None, None, None, None, None, None, None,
         None, "Rembi", 0.88, date.today()),
        
        ("DZ-2024-003", "Rougi", "Hamra", "Femelle", 2.5, "2021-11-10",
         67.8, 78.9, 88.4, 28.7, 61.2, 3,
         22.8, 16.3, 9.8, 1.9, 7, 980, 6,
         None, "Hamra", 0.85, date.today()),
        
        ("DZ-2024-004", "Laitier", "D'man", "Mâle", 5.0, "2019-08-05",
         64.5, 76.3, 84.9, 27.8, 56.8, 4,
         None, None, None, None, None, None, None,
         None, "D'man", 0.90, date.today()),
        
        ("DZ-2024-005", "Kabyle", "Berbère", "Femelle", 3.0, "2021-01-30",
         71.8, 80.2, 92.1, 31.5, 67.9, 4,
         24.1, 17.5, 10.2, 2.0, 8, 1120, 7,
         None, "Berbère", 0.87, date.today())
    ]
    
    for animal in animaux_demo:
        db.execute(
            """INSERT OR IGNORE INTO animaux 
               (identifiant, nom, race, sexe, age_annees, date_naissance,
                hauteur_garrot_cm, longueur_corps_cm, tour_poitrine_cm, largeur_bassin_cm, poids_kg, note_condition,
                largeur_mamelle_cm, hauteur_mamelle_cm, distance_tetines_cm, diametre_tetine_cm, symetrie_score, volume_mamelle_ml, score_morpho_mamelle,
                image_path, race_detectee, confiance_detection, date_analyse) 
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            animal
        )
    
    # Ajouter des données de production laitière
    production_demo = [
        ("DZ-2024-001", date.today() - pd.Timedelta(days=30), 1.8, 6.2, 5.1, 450000),
        ("DZ-2024-001", date.today() - pd.Timedelta(days=15), 2.1, 6.5, 5.3, 380000),
        ("DZ-2024-003", date.today() - pd.Timedelta(days=25), 1.5, 5.8, 4.9, 520000),
        ("DZ-2024-005", date.today() - pd.Timedelta(days=20), 1.9, 6.3, 5.2, 410000)
    ]
    
    for prod in production_demo:
        db.execute(
            """INSERT OR IGNORE INTO production_laitiere 
               (animal_id, date_controle, quantite_lait_l, taux_butyreux, taux_proteique, cellules_somatiques) 
               VALUES (?, ?, ?, ?, ?, ?)""",
            prod
        )

# ============================================================================
# 3. SYSTEME DE MESURE MORPHOMETRIQUE
# ============================================================================

class SystemeMesure:
    """Système de mesure morphométrique avec calibration"""
    
    OBJETS_REFERENCE = {
        "baton_1m": {"longueur_cm": 100.0, "description": "Bâton de 1 mètre"},
        "feuille_a4": {"longueur_cm": 29.7, "description": "Feuille A4 (hauteur)"},
        "carte_bancaire": {"longueur_cm": 8.56, "description": "Carte bancaire standard"},
        "regle_30cm": {"longueur_cm": 30.0, "description": "Règle de 30 cm"},
        "piece_10da": {"diametre_cm": 2.7, "description": "Pièce de 10 DA"}
    }
    
    @staticmethod
    def calculer_echelle(objet_reference, pixels_mesures):
        """Calcule l'échelle pixels → cm"""
        if objet_reference in SystemeMesure.OBJETS_REFERENCE:
            longueur_reelle = SystemeMesure.OBJETS_REFERENCE[objet_reference]["longueur_cm"]
            if pixels_mesures > 0:
                return longueur_reelle / pixels_mesures
        return 0.1  # Valeur par défaut
    
    @staticmethod
    def estimer_poids(tour_poitrine_cm, longueur_corps_cm):
        """Estime le poids à partir des mesures (formule pour ovins)"""
        # Formule simplifiée: Poids (kg) = (Tour de poitrine² × Longueur) / 10800
        return (tour_poitrine_cm ** 2 * longueur_corps_cm) / 10800
    
    @staticmethod
    def calculer_indice_corporel(poids_kg, hauteur_cm):
        """Calcule l'indice corporel"""
        if hauteur_cm > 0:
            return poids_kg / ((hauteur_cm / 100) ** 2)
        return 0
    
    @staticmethod
    def evaluer_mamelle(largeur_cm, hauteur_cm, distance_tetines_cm=None):
        """Évalue la qualité morphologique de la mamelle"""
        score = 5  # Moyen par défaut
        
        # Ratio largeur/hauteur (idéal: 1.3-1.7)
        if hauteur_cm > 0:
            ratio = largeur_cm / hauteur_cm
            if 1.3 <= ratio <= 1.7:
                score += 2
            elif 1.1 <= ratio <= 1.9:
                score += 1
        
        # Distance entre tétines (idéal: 8-12 cm)
        if distance_tetines_cm:
            if 8 <= distance_tetines_cm <= 12:
                score += 1
            elif distance_tetines_cm < 6 or distance_tetines_cm > 15:
                score -= 1
        
        # Limiter entre 1 et 9
        return min(9, max(1, score))
    
    @staticmethod
    def estimer_volume_mamelle(largeur_cm, hauteur_cm):
        """Estime le volume mammaire en ml"""
        # Estimation simplifiée: Volume ≈ largeur × hauteur × 0.6 × 10
        return largeur_cm * hauteur_cm * 0.6 * 10
    
    @staticmethod
    def recommander_race(caracteristiques):
        """Recommande une race basée sur les caractéristiques"""
        recommandations = []
        
        if caracteristiques.get("hauteur", 70) >= 73:
            recommandations.extend(["Ouled Djellal", "Berbère", "Sidaho"])
        elif caracteristiques.get("hauteur", 70) <= 67:
            recommandations.extend(["D'man", "Touareg", "Hamra"])
        else:
            recommandations.extend(["Rembi", "Ouled Djellal", "Berbère"])
        
        return list(set(recommandations))[:3]  # Retourne les 3 premières uniques

# ============================================================================
# 4. MODULE D'ENREGISTREMENT AVEC MESURES
# ============================================================================

def module_enregistrement(db):
    """Module d'enregistrement avec mesures morphométriques"""
    st.header("📝 Enregistrement avec Mesures Morphométriques")
    
    with st.form("form_enregistrement_complet"):
        st.subheader("1. Informations de base")
        
        col_id, col_nom = st.columns(2)
        with col_id:
            identifiant = st.text_input(
                "Identifiant unique*",
                value=f"DZ-{datetime.now().year}-{random.randint(100, 999)}",
                help="Ex: DZ-2024-123"
            )
        with col_nom:
            nom = st.text_input("Nom (optionnel)")
        
        col_race, col_sexe = st.columns(2)
        with col_race:
            race = st.selectbox(
                "Race*",
                db.fetch("SELECT nom FROM races_algeriennes ORDER BY nom")['nom'].tolist()
            )
        with col_sexe:
            sexe = st.selectbox("Sexe*", ["Femelle", "Mâle"])
        
        col_age, col_date = st.columns(2)
        with col_age:
            age = st.number_input("Âge (années)*", min_value=0.0, max_value=15.0, value=3.0, step=0.5)
        with col_date:
            date_naissance = st.date_input("Date de naissance (approximative)")
        
        st.divider()
        st.subheader("2. 📏 Mesures Corporelles")
        
        st.info("""
        **Instructions:**
        1. Utilisez un ruban métrique ou une toise
        2. Mesurez l'animal debout sur surface plane
        3. Pour le tour de poitrine: mesurer derrière les membres antérieurs
        """)
        
        col_mes1, col_mes2 = st.columns(2)
        
        with col_mes1:
            hauteur = st.number_input("Hauteur au garrot (cm)*", 
                                    min_value=30.0, max_value=120.0, value=75.0, step=0.1)
            longueur = st.number_input("Longueur du corps (cm)", 
                                     min_value=30.0, max_value=150.0, value=80.0, step=0.1)
            tour_poitrine = st.number_input("Tour de poitrine (cm)*", 
                                          min_value=40.0, max_value=150.0, value=95.0, step=0.1)
        
        with col_mes2:
            largeur_bassin = st.number_input("Largeur du bassin (cm)", 
                                           min_value=15.0, max_value=60.0, value=30.0, step=0.1)
            poids = st.number_input("Poids (kg) - si connu", 
                                  min_value=10.0, max_value=150.0, value=0.0, step=0.1)
            note_condition = st.slider("Note d'état corporel (1-5)", 1, 5, 3,
                                     help="1: Très maigre, 3: Optimal, 5: Obèse")
        
        # Estimation du poids si non renseigné
        if poids == 0 and tour_poitrine > 0 and longueur > 0:
            poids_estime = SystemeMesure.estimer_poids(tour_poitrine, longueur)
            st.info(f"📊 **Poids estimé:** {poids_estime:.1f} kg (basé sur tour de poitrine et longueur)")
            poids = poids_estime
        
        # Calcul de l'indice corporel
        if poids > 0 and hauteur > 0:
            indice_corporel = SystemeMesure.calculer_indice_corporel(poids, hauteur)
            st.metric("Indice corporel estimé", f"{indice_corporel:.1f}")
        
        st.divider()
        
        # Mesures mammaires pour les femelles
        if sexe == "Femelle":
            st.subheader("3. 🥛 Mesures Mammaires")
            
            col_mam1, col_mam2 = st.columns(2)
            
            with col_mam1:
                largeur_mamelle = st.number_input("Largeur de la mamelle (cm)", 
                                                min_value=10.0, max_value=50.0, value=25.0, step=0.1)
                hauteur_mamelle = st.number_input("Hauteur de la mamelle (cm)", 
                                                min_value=5.0, max_value=40.0, value=18.0, step=0.1)
            
            with col_mam2:
                distance_tetines = st.number_input("Distance entre tétines (cm)", 
                                                 min_value=5.0, max_value=30.0, value=10.5, step=0.1)
                diametre_tetine = st.number_input("Diamètre des tétines (cm)", 
                                                min_value=0.5, max_value=5.0, value=2.0, step=0.1)
            
            # Calculs automatiques
            if largeur_mamelle > 0 and hauteur_mamelle > 0:
                score_mamelle = SystemeMesure.evaluer_mamelle(largeur_mamelle, hauteur_mamelle, distance_tetines)
                volume_mamelle = SystemeMesure.estimer_volume_mamelle(largeur_mamelle, hauteur_mamelle)
                
                col_score, col_volume = st.columns(2)
                with col_score:
                    # Afficher le score avec couleur
                    if score_mamelle >= 7:
                        couleur = "🟢"
                        appreciation = "Excellente"
                    elif score_mamelle >= 5:
                        couleur = "🟡"
                        appreciation = "Bonne"
                    else:
                        couleur = "🔴"
                        appreciation = "À améliorer"
                    
                    st.metric("Score morphologique mamelle", f"{couleur} {score_mamelle}/9", appreciation)
                
                with col_volume:
                    st.metric("Volume mammaire estimé", f"{volume_mamelle:.0f} ml")
        
        st.divider()
        st.subheader("4. Calibration avec référence (optionnel)")
        
        avec_calibration = st.checkbox("Utiliser un objet de référence pour conversion photo")
        
        if avec_calibration:
            col_ref1, col_ref2 = st.columns(2)
            
            with col_ref1:
                objet_ref = st.selectbox(
                    "Objet de référence",
                    list(SystemeMesure.OBJETS_REFERENCE.keys())
                )
                ref_info = SystemeMesure.OBJETS_REFERENCE[objet_ref]
                st.info(f"**{ref_info['description']}** ({ref_info['longueur_cm']} cm)")
            
            with col_ref2:
                pixels_reference = st.number_input(
                    "Taille de la référence sur la photo (pixels)",
                    min_value=1.0,
                    max_value=10000.0,
                    value=500.0,
                    step=10.0,
                    help="Mesurez l'objet de référence sur la photo avec un logiciel"
                )
                
                if pixels_reference > 0:
                    echelle = ref_info["longueur_cm"] / pixels_reference
                    st.success(f"Échelle: **1 pixel = {echelle:.4f} cm**")
        
        # Bouton de soumission
        submitted = st.form_submit_button("✅ Enregistrer l'animal avec mesures", type="primary")
        
        if submitted:
            # Validation
            if not identifiant or not race:
                st.error("Veuillez remplir tous les champs obligatoires (*)")
                return
            
            try:
                # Préparer les données
                animal_data = {
                    'identifiant': identifiant,
                    'nom': nom,
                    'race': race,
                    'sexe': sexe,
                    'age_annees': age,
                    'date_naissance': date_naissance,
                    'hauteur_garrot_cm': hauteur,
                    'longueur_corps_cm': longueur,
                    'tour_poitrine_cm': tour_poitrine,
                    'largeur_bassin_cm': largeur_bassin,
                    'poids_kg': poids,
                    'note_condition': note_condition,
                    'race_detectee': race,
                    'confiance_detection': 1.0,
                    'date_analyse': date.today()
                }
                
                # Ajouter les mesures mammaires pour femelles
                if sexe == "Femelle":
                    score_mamelle = SystemeMesure.evaluer_mamelle(largeur_mamelle, hauteur_mamelle, distance_tetines)
                    volume_mamelle = SystemeMesure.estimer_volume_mamelle(largeur_mamelle, hauteur_mamelle)
                    
                    animal_data.update({
                        'largeur_mamelle_cm': largeur_mamelle,
                        'hauteur_mamelle_cm': hauteur_mamelle,
                        'distance_tetines_cm': distance_tetines,
                        'diametre_tetine_cm': diametre_tetine,
                        'symetrie_score': 8,  # Valeur par défaut
                        'volume_mamelle_ml': volume_mamelle,
                        'score_morpho_mamelle': score_mamelle
                    })
                
                # Construire la requête SQL
                columns = ', '.join(animal_data.keys())
                placeholders = ', '.join(['?' for _ in animal_data])
                values = tuple(animal_data.values())
                
                query = f"INSERT OR REPLACE INTO animaux ({columns}) VALUES ({placeholders})"
                
                # Exécution
                db.execute(query, values)
                
                st.success(f"✅ Animal {identifiant} enregistré avec succès !")
                st.balloons()
                
                # Afficher un résumé
                with st.expander("📋 Voir le résumé des mesures"):
                    st.json(animal_data)
                    
                # Recommandations basées sur les mesures
                st.subheader("💡 Recommandations basées sur les mesures")
                
                recommandations = SystemeMesure.recommander_race({
                    "hauteur": hauteur,
                    "poids": poids,
                    "tour_poitrine": tour_poitrine
                })
                
                if race in recommandations:
                    st.success(f"✅ La race {race} est bien adaptée aux mesures de cet animal.")
                else:
                    st.warning(f"⚠️ La race {race} pourrait ne pas être optimale pour ces mesures.")
                    st.info(f"Races recommandées: {', '.join(recommandations)}")
                
            except Exception as e:
                st.error(f"❌ Erreur lors de l'enregistrement: {str(e)}")

# ============================================================================
# 5. MODULE ANALYSE ET STATISTIQUES
# ============================================================================

def module_analyses(db):
    """Module d'analyses et statistiques"""
    st.header("📊 Analyses et Statistiques du Troupeau")
    
    # Récupérer toutes les données
    animaux = db.fetch("""
        SELECT * FROM animaux 
        ORDER BY date_enregistrement DESC
    """)
    
    if animaux.empty:
        st.info("Aucun animal enregistré. Ajoutez des animaux pour voir les analyses.")
        return
    
    # KPI Principaux
    st.subheader("📈 Indicateurs Clés")
    
    col_kpi1, col_kpi2, col_kpi3, col_kpi4 = st.columns(4)
    
    with col_kpi1:
        total = len(animaux)
        st.metric("Total animaux", total)
    
    with col_kpi2:
        femelles = len(animaux[animaux['sexe'] == 'Femelle'])
        st.metric("Femelles", femelles)
    
    with col_kpi3:
        if 'poids_kg' in animaux.columns:
            poids_moyen = animaux['poids_kg'].mean()
            st.metric("Poids moyen", f"{poids_moyen:.1f} kg")
    
    with col_kpi4:
        if 'hauteur_garrot_cm' in animaux.columns:
            hauteur_moyenne = animaux['hauteur_garrot_cm'].mean()
            st.metric("Hauteur moyenne", f"{hauteur_moyenne:.1f} cm")
    
    # Analyse par race
    st.subheader("🏷️ Analyse par Race")
    
    if 'race' in animaux.columns and not animaux.empty:
        # Statistiques par race
        stats_race = animaux.groupby('race').agg({
            'hauteur_garrot_cm': ['mean', 'std', 'count'],
            'poids_kg': ['mean', 'std'],
            'tour_poitrine_cm': 'mean'
        }).round(1)
        
        # Formater le DataFrame
        stats_race.columns = ['_'.join(col).strip() for col in stats_race.columns.values]
        stats_race = stats_race.rename(columns={
            'hauteur_garrot_cm_mean': 'Hauteur moyenne (cm)',
            'hauteur_garrot_cm_std': 'Écart-type hauteur',
            'hauteur_garrot_cm_count': 'Nombre',
            'poids_kg_mean': 'Poids moyen (kg)',
            'poids_kg_std': 'Écart-type poids',
            'tour_poitrine_cm_mean': 'Tour poitrine moyen (cm)'
        })
        
        st.dataframe(stats_race)
    
    # Graphiques
    st.subheader("📈 Visualisations")
    
    col_graph1, col_graph2 = st.columns(2)
    
    with col_graph1:
        if 'race' in animaux.columns:
            import plotly.express as px
            
            # Distribution des races
            fig1 = px.pie(
                animaux, 
                names='race',
                title='Répartition par race',
                hole=0.3,
                color_discrete_sequence=px.colors.sequential.RdBu
            )
            st.plotly_chart(fig1, use_container_width=True)
    
    with col_graph2:
        if all(col in animaux.columns for col in ['hauteur_garrot_cm', 'poids_kg']):
            import plotly.express as px
            
            # Relation hauteur-poids
            fig2 = px.scatter(
                animaux,
                x='hauteur_garrot_cm',
                y='poids_kg',
                color='race',
                size='tour_poitrine_cm',
                title='Relation Hauteur-Poids par Race',
                hover_data=['identifiant', 'nom', 'age_annees'],
                trendline="ols"
            )
            st.plotly_chart(fig2, use_container_width=True)
    
    # Analyse mammaire (si femelles)
    st.subheader("🥛 Analyse des Mamelles (Femelles)")
    
    femelles = animaux[animaux['sexe'] == 'Femelle']
    if not femelles.empty and 'score_morpho_mamelle' in femelles.columns:
        col_mam1, col_mam2, col_mam3 = st.columns(3)
        
        with col_mam1:
            score_moyen = femelles['score_morpho_mamelle'].mean()
            st.metric("Score mammaire moyen", f"{score_moyen:.1f}/9")
        
        with col_mam2:
            if 'largeur_mamelle_cm' in femelles.columns:
                largeur_moyenne = femelles['largeur_mamelle_cm'].mean()
                st.metric("Largeur mamelle moyenne", f"{largeur_moyenne:.1f} cm")
        
        with col_mam3:
            if 'volume_mamelle_ml' in femelles.columns:
                volume_moyen = femelles['volume_mamelle_ml'].mean()
                st.metric("Volume mammaire moyen", f"{volume_moyen:.0f} ml")
        
        # Top 5 des meilleures mamelles
        st.subheader("🏆 Top 5 des meilleures mamelles")
        top_mamelles = femelles.nlargest(5, 'score_morpho_mamelle')[['identifiant', 'nom', 'race', 'score_morpho_mamelle', 'largeur_mamelle_cm', 'volume_mamelle_ml']]
        st.dataframe(top_mamelles)
    
    # Export des données
    st.subheader("💾 Export des données")
    
    col_exp1, col_exp2 = st.columns(2)
    
    with col_exp1:
        if st.button("📥 Exporter en CSV"):
            csv = animaux.to_csv(index=False)
            st.download_button(
                label="Télécharger CSV",
                data=csv,
                file_name=f"troupeau_ovin_{date.today()}.csv",
                mime="text/csv"
            )
    
    with col_exp2:
        if st.button("📊 Générer rapport complet"):
            # Créer un rapport simplifié
            rapport = f"""
            RAPPORT DU TROUPEAU OVIN - {date.today()}
            =========================================
            
            Nombre total d'animaux: {len(animaux)}
            Nombre de femelles: {len(femelles)}
            Nombre de mâles: {len(animaux) - len(femelles)}
            
            """
            
            if 'poids_kg' in animaux.columns:
                rapport += f"Poids moyen: {animaux['poids_kg'].mean():.1f} kg\n"
            
            if 'hauteur_garrot_cm' in animaux.columns:
                rapport += f"Hauteur moyenne: {animaux['hauteur_garrot_cm'].mean():.1f} cm\n"
            
            st.download_button(
                label="Télécharger rapport",
                data=rapport,
                file_name=f"rapport_troupeau_{date.today()}.txt",
                mime="text/plain"
            )

# ============================================================================
# 6. MODULE RECHERCHE AVANCEE
# ============================================================================

def module_recherche(db):
    """Module de recherche avancée"""
    st.header("🔍 Recherche Avancée")
    
    # Filtres
    st.subheader("Filtres de recherche")
    
    col_filtre1, col_filtre2 = st.columns(2)
    
    with col_filtre1:
        race_filtre = st.selectbox(
            "Race",
            ["Toutes"] + db.fetch("SELECT DISTINCT nom FROM races_algeriennes ORDER BY nom")['nom'].tolist()
        )
        
        sexe_filtre = st.selectbox(
            "Sexe",
            ["Tous", "Femelle", "Mâle"]
        )
    
    with col_filtre2:
        age_min = st.slider("Âge minimum (ans)", 0, 20, 0)
        age_max = st.slider("Âge maximum (ans)", 0, 20, 20)
        
        if st.checkbox("Filtrer par mesures"):
            hauteur_min = st.number_input("Hauteur minimum (cm)", 30.0, 120.0, 60.0)
            hauteur_max = st.number_input("Hauteur maximum (cm)", 30.0, 120.0, 85.0)
    
    # Construire la requête
    query = "SELECT * FROM animaux WHERE 1=1"
    params = []
    
    if race_filtre != "Toutes":
        query += " AND race = ?"
        params.append(race_filtre)
    
    if sexe_filtre != "Tous":
        query += " AND sexe = ?"
        params.append(sexe_filtre)
    
    if age_min > 0:
        query += " AND age_annees >= ?"
        params.append(age_min)
    
    if age_max < 20:
        query += " AND age_annees <= ?"
        params.append(age_max)
    
    # Exécuter la recherche
    resultats = db.fetch(query, tuple(params) if params else ())
    
    if not resultats.empty:
        st.success(f"✅ {len(resultats)} animaux trouvés")
        
        # Options d'affichage
        affichage = st.radio(
            "Mode d'affichage",
            ["Tableau complet", "Vue résumée", "Cartes individuelles"]
        )
        
        if affichage == "Tableau complet":
            st.dataframe(resultats)
        
        elif affichage == "Vue résumée":
            colonnes_resume = ['identifiant', 'nom', 'race', 'sexe', 'age_annees', 
                             'hauteur_garrot_cm', 'poids_kg', 'note_condition']
            colonnes_dispo = [col for col in colonnes_resume if col in resultats.columns]
            st.dataframe(resultats[colonnes_dispo])
        
        else:  # Cartes individuelles
            for _, animal in resultats.iterrows():
                with st.expander(f"{animal['identifiant']} - {animal['nom'] or 'Sans nom'}"):
                    col_card1, col_card2 = st.columns(2)
                    
                    with col_card1:
                        st.metric("Race", animal['race'])
                        st.metric("Sexe", animal['sexe'])
                        st.metric("Âge", f"{animal['age_annees']} ans")
                    
                    with col_card2:
                        if pd.notna(animal['hauteur_garrot_cm']):
                            st.metric("Hauteur", f"{animal['hauteur_garrot_cm']} cm")
                        if pd.notna(animal['poids_kg']):
                            st.metric("Poids", f"{animal['poids_kg']} kg")
                        if pd.notna(animal['note_condition']):
                            st.metric("État corporel", f"{animal['note_condition']}/5")
        
        # Statistiques des résultats
        st.subheader("📊 Statistiques des résultats")
        
        if len(resultats) > 1:
            col_stats1, col_stats2 = st.columns(2)
            
            with col_stats1:
                if 'hauteur_garrot_cm' in resultats.columns:
                    avg_height = resultats['hauteur_garrot_cm'].mean()
                    st.metric("Hauteur moyenne", f"{avg_height:.1f} cm")
                
                if 'age_annees' in resultats.columns:
                    avg_age = resultats['age_annees'].mean()
                    st.metric("Âge moyen", f"{avg_age:.1f} ans")
            
            with col_stats2:
                if 'poids_kg' in resultats.columns:
                    avg_weight = resultats['poids_kg'].mean()
                    st.metric("Poids moyen", f"{avg_weight:.1f} kg")
    else:
        st.info("Aucun animal ne correspond aux critères de recherche.")

# ============================================================================
# 7. INTERFACE PRINCIPALE
# ============================================================================

def main():
    """Fonction principale de l'application"""
    
    # Initialisation de la base de données
    db = DatabaseManager()
    init_database(db)
    
    # Sidebar
    with st.sidebar:
        st.markdown("# 🐑 Expert Ovin DZ")
        st.markdown("---")
        
        # Menu de navigation
        menu = st.radio(
            "Navigation",
            [
                "🏠 Tableau de bord",
                "📝 Enregistrement avec mesures",
                "🔍 Recherche avancée",
                "📊 Analyses et statistiques",
                "⚙️ Configuration"
            ]
        )
        
        st.markdown("---")
        
        # Boutons d'action
        col_btn1, col_btn2 = st.columns(2)
        
        with col_btn1:
            if st.button("🔄 Démo", use_container_width=True):
                with st.spinner("Création des données de démo..."):
                    creer_donnees_demo(db)
                    st.success("✅ Données de démo créées !")
        
        with col_btn2:
            if st.button("🗑️ Nettoyer", use_container_width=True):
                if st.checkbox("Confirmer la suppression des données"):
                    db.execute("DELETE FROM animaux")
                    db.execute("DELETE FROM production_laitiere")
                    db.execute("DELETE FROM gestations")
                    st.warning("Base de données nettoyée !")
        
        st.markdown("---")
        st.caption(f"Version 1.0 | {date.today()}")
    
    # Contenu principal selon le menu
    if menu == "🏠 Tableau de bord":
        st.header("🏠 Tableau de bord")
        
        # Statistiques rapides
        col_dash1, col_dash2, col_dash3, col_dash4 = st.columns(4)
        
        with col_dash1:
            total = db.fetch("SELECT COUNT(*) as count FROM animaux")['count'][0]
            st.metric("Animaux total", total)
        
        with col_dash2:
            races = db.fetch("SELECT COUNT(DISTINCT race) as count FROM animaux")['count'][0]
            st.metric("Races", races)
        
        with col_dash3:
            femelles = db.fetch("SELECT COUNT(*) as count FROM animaux WHERE sexe = 'Femelle'")['count'][0]
            st.metric("Femelles", femelles)
        
        with col_dash4:
            mâles = total - femelles
            st.metric("Mâles", mâles)
        
        # Derniers enregistrements
        st.subheader("🆕 Derniers animaux enregistrés")
        
        derniers = db.fetch("""
            SELECT identifiant, nom, race, sexe, age_annees, 
                   hauteur_garrot_cm, poids_kg, date_enregistrement
            FROM animaux
            ORDER BY date_enregistrement DESC
            LIMIT 5
        """)
        
        if not derniers.empty:
            st.dataframe(derniers)
        else:
            st.info("Aucun animal enregistré. Commencez par en ajouter !")
        
        # Guide de démarrage
        st.subheader("🚀 Guide de démarrage rapide")
        
        col_guide1, col_guide2 = st.columns(2)
        
        with col_guide1:
            st.info("""
            **Pour commencer:**
            1. Enregistrez vos animaux
            2. Prenez leurs mesures
            3. Consultez les analyses
            """)
        
        with col_guide2:
            st.info("""
            **Mesures importantes:**
            - Hauteur au garrot
            - Tour de poitrine
            - Poids ou estimation
            - Mesures mammaires (femelles)
            """)
    
    elif menu == "📝 Enregistrement avec mesures":
        module_enregistrement(db)
    
    elif menu == "🔍 Recherche avancée":
        module_recherche(db)
    
    elif menu == "📊 Analyses et statistiques":
        module_analyses(db)
    
    elif menu == "⚙️ Configuration":
        st.header("⚙️ Configuration")
        
        st.subheader("Base de données")
        
        if st.button("🔍 Vérifier l'état de la base"):
            tables = db.fetch("""
                SELECT name, sql 
                FROM sqlite_master 
                WHERE type='table'
            """)
            
            if not tables.empty:
                st.success(f"✅ Base de données opérationnelle ({len(tables)} tables)")
                
                for _, row in tables.iterrows():
                    with st.expander(f"Table: {row['name']}"):
                        st.code(row['sql'])
            else:
                st.error("❌ Base de données vide ou corrompue")
        
        st.subheader("Installation requise")
        
        st.code("""
# Packages nécessaires pour cette version:
streamlit==1.28.0
pandas==2.1.4
numpy==1.24.3
plotly==5.18.0
openpyxl==3.1.2
Pillow==10.1.0
python-dateutil==2.8.2
requests==2.31.0
        """)
        
        st.info("""
        **Configuration minimale:**
        - Python 3.8+
        - 1GB RAM
        - 50MB espace disque
        
        **Fonctionnalités incluses:**
        - Base de données SQLite
        - Mesures morphométriques
        - Analyse statistique
        - Export des données
        - Interface graphique
        """)

# ============================================================================
# POINT D'ENTREE
# ============================================================================

if __name__ == "__main__":
    # Vérification des imports critiques
    try:
        import streamlit
        import pandas
        import numpy
        st.success("✅ Packages principaux chargés avec succès")
        main()
    except ImportError as e:
        st.error(f"❌ Package manquant: {e}")
        st.info("Veuillez installer les packages requis avec le fichier requirements.txt fourni")
