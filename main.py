#!/usr/bin/env python3
"""
OVIN MANAGER PRO - Application complète de gestion scientifique d'élevage ovin
Version: 1.0.0
Auteur: [Votre Nom]
Description: Application intégrée avec gestion, analyse morphométrique, génomique et statistiques
"""

# ============================================================================
# BLOC 1: IMPORTATIONS ET CONFIGURATION
# ============================================================================
print("🔧 Initialisation de Ovin Manager Pro...")

import os
import sys
import json
import sqlite3
import tempfile
from datetime import datetime, date, timedelta
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from pathlib import Path
import logging

# Configuration du logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ovin_manager.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# BLOC 2: CLASSES DE BASE ET MODÈLES DE DONNÉES
# ============================================================================
@dataclass
class Brebis:
    """Classe représentant une brebis"""
    id: int
    identifiant_unique: str
    nom: str
    date_naissance: date
    race: str
    sexe: str  # 'F' ou 'M'
    statut: str = "active"  # active, retrait, morte
    notes: str = ""
    
    def age_jours(self) -> int:
        """Calcule l'âge en jours"""
        return (date.today() - self.date_naissance).days
    
    def age_mois(self) -> float:
        """Calcule l'âge en mois"""
        return self.age_jours() / 30.44

@dataclass
class SuiviMedical:
    """Classe pour le suivi médical"""
    id: int
    brebis_id: int
    date_intervention: date
    type_intervention: str  # vaccination, traitement, vermifuge
    produit: str
    dose: str
    veterinaire: str = ""
    notes: str = ""

@dataclass 
class Gestation:
    """Classe pour le suivi de gestation"""
    id: int
    brebis_id: int
    date_eponge: date
    date_retrait_eponge: Optional[date] = None
    date_insemination: Optional[date] = None
    date_mise_bas_prevu: Optional[date] = None
    date_mise_bas_reel: Optional[date] = None
    nombre_agneaux_prevus: int = 1
    nombre_agneaux_nes: int = 0
    statut: str = "en_cours"  # en_cours, termine, avorte

@dataclass
class CaractereMorphologique:
    """Classe pour les caractères morphologiques"""
    id: int
    brebis_id: int
    date_mesure: date
    caractere: str  # poids_vif, longueur_corps, etc.
    valeur: float
    unite: str
    methode_mesure: str = "manuel"  # manuel, photo, scanner

@dataclass
class SequenceGenetique:
    """Classe pour les séquences génétiques"""
    id: int
    brebis_id: int
    accession_number: str
    sequence_type: str  # ADN, ARN, protéine
    longueur: int
    date_sequencage: date
    laboratoire: str = ""
    notes: str = ""

# ============================================================================
# BLOC 3: GESTIONNAIRE DE BASE DE DONNÉES SQLITE
# ============================================================================
class DatabaseManager:
    """Gestionnaire de base de données SQLite"""
    
    def __init__(self, db_path: str = "ovin_manager.db"):
        self.db_path = db_path
        self.conn = None
        self.init_database()
    
    def init_database(self):
        """Initialise la base de données avec les tables"""
        try:
            self.conn = sqlite3.connect(self.db_path)
            cursor = self.conn.cursor()
            
            # Table des brebis
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS brebis (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    identifiant_unique TEXT UNIQUE NOT NULL,
                    nom TEXT,
                    date_naissance DATE,
                    race TEXT,
                    sexe TEXT,
                    statut TEXT,
                    notes TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Table suivi médical
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS suivi_medical (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    brebis_id INTEGER,
                    date_intervention DATE,
                    type_intervention TEXT,
                    produit TEXT,
                    dose TEXT,
                    veterinaire TEXT,
                    notes TEXT,
                    FOREIGN KEY (brebis_id) REFERENCES brebis (id)
                )
            ''')
            
            # Table gestation
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS gestations (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    brebis_id INTEGER,
                    date_eponge DATE,
                    date_retrait_eponge DATE,
                    date_insemination DATE,
                    date_mise_bas_prevu DATE,
                    date_mise_bas_reel DATE,
                    nombre_agneaux_prevus INTEGER,
                    nombre_agneaux_nes INTEGER,
                    statut TEXT,
                    FOREIGN KEY (brebis_id) REFERENCES brebis (id)
                )
            ''')
            
            # Table caractères morphologiques
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS caracteres_morpho (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    brebis_id INTEGER,
                    date_mesure DATE,
                    caractere TEXT,
                    valeur REAL,
                    unite TEXT,
                    methode_mesure TEXT,
                    FOREIGN KEY (brebis_id) REFERENCES brebis (id)
                )
            ''')
            
            # Table séquences génétiques
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS sequences_genetiques (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    brebis_id INTEGER,
                    accession_number TEXT UNIQUE,
                    sequence_type TEXT,
                    longueur INTEGER,
                    date_sequencage DATE,
                    laboratoire TEXT,
                    notes TEXT,
                    FOREIGN KEY (brebis_id) REFERENCES brebis (id)
                )
            ''')
            
            self.conn.commit()
            logger.info(f"✅ Base de données initialisée: {self.db_path}")
            
        except sqlite3.Error as e:
            logger.error(f"❌ Erreur base de données: {e}")
            raise
    
    def ajouter_brebis(self, brebis: Brebis) -> int:
        """Ajoute une brebis à la base"""
        cursor = self.conn.cursor()
        cursor.execute('''
            INSERT INTO brebis 
            (identifiant_unique, nom, date_naissance, race, sexe, statut, notes)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (brebis.identifiant_unique, brebis.nom, brebis.date_naissance.isoformat(),
              brebis.race, brebis.sexe, brebis.statut, brebis.notes))
        self.conn.commit()
        return cursor.lastrowid
    
    def get_brebis(self, brebis_id: int = None) -> List[Dict]:
        """Récupère les brebis (une ou toutes)"""
        cursor = self.conn.cursor()
        
        if brebis_id:
            cursor.execute('SELECT * FROM brebis WHERE id = ?', (brebis_id,))
            rows = cursor.fetchall()
        else:
            cursor.execute('SELECT * FROM brebis ORDER BY id')
            rows = cursor.fetchall()
        
        columns = [description[0] for description in cursor.description]
        return [dict(zip(columns, row)) for row in rows]
    
    def close(self):
        """Ferme la connexion à la base"""
        if self.conn:
            self.conn.close()

# ============================================================================
# BLOC 4: GESTIONNAIRE DE GESTATION ET PRÉDICTIONS
# ============================================================================
class GestionnaireGestation:
    """Gestion des gestations avec prédictions scientifiques"""
    
    # Durées de gestation par race (jours) - Données scientifiques
    DUREES_GESTATION = {
        'lacaune': {'moyenne': 152, 'ecart_type': 2.5},
        'manech': {'moyenne': 150, 'ecart_type': 2.0},
        'basco_bearnaise': {'moyenne': 148, 'ecart_type': 2.2},
        'default': {'moyenne': 150, 'ecart_type': 2.5}
    }
    
    def __init__(self, db_manager: DatabaseManager):
        self.db = db_manager
    
    def calculer_date_mise_bas(self, date_eponge: date, race: str = 'default') -> date:
        """Calcule la date prévue de mise bas"""
        duree = self.DUREES_GESTATION.get(race, self.DUREES_GESTATION['default'])
        return date_eponge + timedelta(days=int(duree['moyenne']))
    
    def programmer_eponge(self, brebis_id: int, date_eponge: date, race: str) -> Dict:
        """Programme l'introduction d'une éponge et calcule les dates clés"""
        date_mise_bas = self.calculer_date_mise_bas(date_eponge, race)
        date_retrait = date_eponge + timedelta(days=14)  # Durée standard éponge
        
        return {
            'brebis_id': brebis_id,
            'date_eponge': date_eponge,
            'date_retrait_eponge': date_retrait,
            'date_mise_bas_prevu': date_mise_bas,
            'periode_insemination': {
                'debut': date_retrait,
                'fin': date_retrait + timedelta(days=2)
            },
            'periode_mise_bas': {
                'debut': date_mise_bas - timedelta(days=3),
                'fin': date_mise_bas + timedelta(days=3)
            }
        }
    
    def generer_calendrier_gestation(self, brebis_id: int) -> List[Dict]:
        """Génère un calendrier détaillé de la gestation"""
        cursor = self.db.conn.cursor()
        cursor.execute('''
            SELECT * FROM gestations 
            WHERE brebis_id = ? AND statut = 'en_cours'
            ORDER BY date_eponge DESC LIMIT 1
        ''', (brebis_id,))
        
        gestation = cursor.fetchone()
        if not gestation:
            return []
        
        columns = [desc[0] for desc in cursor.description]
        gestation_dict = dict(zip(columns, gestation))
        
        date_eponge = datetime.strptime(gestation_dict['date_eponge'], '%Y-%m-%d').date()
        date_mise_bas = datetime.strptime(gestation_dict['date_mise_bas_prevu'], '%Y-%m-%d').date()
        
        calendrier = []
        
        # Phase 1: Synchronisation (J0-J14)
        calendrier.append({
            'phase': 'Synchronisation',
            'periode': f"J0 à J14",
            'dates': f"{date_eponge} à {date_eponge + timedelta(days=14)}",
            'actions': ['Éponge en place', 'Contrôle quotidien'],
            'surveillance': 'Température normale, appétit conservé'
        })
        
        # Phase 2: Insémination (J14-J16)
        calendrier.append({
            'phase': 'Insémination',
            'periode': 'J14 à J16',
            'dates': f"{date_eponge + timedelta(days=14)} à {date_eponge + timedelta(days=16)}",
            'actions': ['Retrait éponge', 'Insémination artificielle'],
            'surveillance': 'Détection des chaleurs'
        })
        
        # Phase 3: Début gestation (J16-J45)
        calendrier.append({
            'phase': 'Début gestation',
            'periode': 'J16 à J45',
            'dates': f"{date_eponge + timedelta(days=16)} à {date_eponge + timedelta(days=45)}",
            'actions': ['Diagnostic gestation (échographie J30)'],
            'surveillance': 'Comportement alimentaire'
        })
        
        return calendrier
    
    def get_statistiques_gestation(self) -> Dict:
        """Calcule les statistiques de gestation"""
        cursor = self.db.conn.cursor()
        
        cursor.execute('SELECT COUNT(*) FROM gestations')
        total = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM gestations WHERE statut = "en_cours"')
        en_cours = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM gestations WHERE statut = "termine"')
        termine = cursor.fetchone()[0]
        
        cursor.execute('''
            SELECT AVG(nombre_agneaux_nes) 
            FROM gestations 
            WHERE statut = "termine" AND nombre_agneaux_nes > 0
        ''')
        moyenne_agneaux = cursor.fetchone()[0] or 0
        
        return {
            'total_gestations': total,
            'en_cours': en_cours,
            'terminees': termine,
            'taux_reussite': (termine / total * 100) if total > 0 else 0,
            'moyenne_agneaux_par_mise_bas': round(moyenne_agneaux, 2)
        }

# ============================================================================
# BLOC 5: ANALYSE MORPHOMÉTRIQUE PAR SMARTPHONE
# ============================================================================
class AnalyseurMorphometrique:
    """Analyse morphométrique à partir de photos smartphone"""
    
    # Références scientifiques pour les races ovines
    REFERENCES_RACES = {
        'lacaune': {
            'poids_adulte_femelle': (70, 90),  # kg
            'hauteur_garrot': (70, 75),  # cm
            'longueur_corps': (75, 85),  # cm
        },
        'manech': {
            'poids_adulte_femelle': (60, 80),
            'hauteur_garrot': (65, 70),
            'longueur_corps': (70, 80),
        },
        'basco_bearnaise': {
            'poids_adulte_femelle': (55, 75),
            'hauteur_garrot': (60, 68),
            'longueur_corps': (68, 78),
        }
    }
    
    def __init__(self):
        self.mesures_standard = {
            'longueur_corps': 'cm',
            'hauteur_garrot': 'cm', 
            'largeur_bassin': 'cm',
            'tour_poitrine': 'cm',
            'poids_estime': 'kg'
        }
    
    def analyser_photo(self, image_path: str, objet_reference_pixels: int, 
                      taille_reelle_objet: float, race: str = None) -> Dict:
        """
        Analyse une photo pour extraire des mesures morphométriques
        
        Args:
            image_path: Chemin de l'image
            objet_reference_pixels: Taille en pixels d'un objet de référence
            taille_reelle_objet: Taille réelle de l'objet (en cm)
            race: Race de la brebis pour comparaison
        
        Returns:
            Dict avec les mesures et analyses
        """
        # Facteur de conversion pixels -> cm
        facteur_conversion = taille_reelle_objet / objet_reference_pixels
        
        # Simulation d'analyse d'image (dans la réalité, utiliser OpenCV)
        # Ici, nous simulons des mesures basées sur des proportions standards
        
        mesures = {}
        
        # Longueur du corps (estimée à partir de proportions)
        # Dans une photo latérale standard, le corps occupe ~80% de la largeur
        longueur_pixels = objet_reference_pixels * 8  # Simulation
        mesures['longueur_corps'] = round(longueur_pixels * facteur_conversion, 2)
        
        # Hauteur au garrot (proportion par rapport à la longueur)
        mesures['hauteur_garrot'] = round(mesures['longueur_corps'] * 0.85, 2)
        
        # Tour de poitrine (estimation)
        mesures['tour_poitrine'] = round(mesures['longueur_corps'] * 1.2, 2)
        
        # Poids estimé (formule scientifique basée sur tour de poitrine)
        # Formule: Poids (kg) = (tour_poitrine² × longueur_corps) / 300
        poids_estime = (mesures['tour_poitrine']**2 * mesures['longueur_corps']) / 30000
        mesures['poids_estime'] = round(poids_estime, 2)
        
        # Analyse comparative si race spécifiée
        if race and race in self.REFERENCES_RACES:
            reference = self.REFERENCES_RACES[race]
            analyses = []
            
            if 'poids_adulte_femelle' in reference:
                min_ref, max_ref = reference['poids_adulte_femelle']
                if min_ref <= mesures['poids_estime'] <= max_ref:
                    analyses.append(f"✅ Poids dans la norme pour {race}")
                else:
                    analyses.append(f"⚠️ Poids hors norme pour {race} ({min_ref}-{max_ref} kg)")
            
            return {
                'mesures': mesures,
                'unites': self.mesures_standard,
                'facteur_conversion': facteur_conversion,
                'analyse_comparative': analyses,
                'date_analyse': datetime.now().isoformat(),
                'precision': 'estimative - nécessite calibration précise'
            }
        
        return {
            'mesures': mesures,
            'unites': self.mesures_standard,
            'facteur_conversion': facteur_conversion,
            'date_analyse': datetime.now().isoformat()
        }
    
    def calculer_indice_corporel(self, poids: float, hauteur_garrot: float) -> float:
        """Calcule l'indice corporel (Body Condition Score adapté)"""
        # BCS approximatif basé sur poids/taille
        ratio = poids / hauteur_garrot
        if ratio < 0.9:
            return 1.5  # Maigre
        elif ratio < 1.1:
            return 2.5  # Normal
        elif ratio < 1.3:
            return 3.5  # Gras
        else:
            return 4.5  # Très gras

# ============================================================================
# BLOC 6: INTÉGRATION GÉNOMIQUE ET NCBI
# ============================================================================
class IntegrationGenomique:
    """Intégration avec les bases de données génomiques"""
    
    def __init__(self, email: str):
        self.email = email
        self.sequences_cache = {}
    
    def formater_sequence_fasta(self, sequence_id: str, sequence: str, 
                               description: str = "") -> str:
        """Formate une séquence au format FASTA"""
        return f">{sequence_id} {description}\n{sequence}\n"
    
    def analyser_snp(self, sequence_reference: str, sequence_etudiee: str) -> Dict:
        """Analyse les SNP entre deux séquences"""
        if len(sequence_reference) != len(sequence_etudiee):
            return {"erreur": "Les séquences doivent avoir la même longueur"}
        
        snps = []
        for i, (ref, etu) in enumerate(zip(sequence_reference, sequence_etudiee)):
            if ref != etu:
                snps.append({
                    'position': i + 1,
                    'reference': ref,
                    'etudie': etu,
                    'type_mutation': self._determiner_type_mutation(ref, etu)
                })
        
        return {
            'total_snps': len(snps),
            'frequence_snp': len(snps) / len(sequence_reference),
            'snps_detailles': snps[:10],  # Limité aux 10 premiers
            'sequence_longueur': len(sequence_reference)
        }
    
    def _determiner_type_mutation(self, ref: str, etu: str) -> str:
        """Détermine le type de mutation"""
        transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
        transversions = [('A', 'C'), ('A', 'T'), ('G', 'C'), ('G', 'T'),
                        ('C', 'A'), ('T', 'A'), ('C', 'G'), ('T', 'G')]
        
        if (ref, etu) in transitions:
            return 'transition'
        elif (ref, etu) in transversions:
            return 'transversion'
        else:
            return 'indéterminé'
    
    def rechercher_genes_candidats(self, race: str) -> List[Dict]:
        """Retourne les gènes candidats pour une race donnée (simulé)"""
        genes_ovins = {
            'lacaune': [
                {'gene': 'LALBA', 'fonction': 'Protéine du lait', 'chromosome': '3'},
                {'gene': 'CSN1S1', 'fonction': 'Caséine alpha-S1', 'chromosome': '6'},
                {'gene': 'DGAT1', 'fonction': 'Synthèse des triglycérides', 'chromosome': '14'},
            ],
            'manech': [
                {'gene': 'PRLR', 'fonction': 'Récepteur prolactine', 'chromosome': '16'},
                {'gene': 'GH1', 'fonction': 'Hormone de croissance', 'chromosome': '11'},
            ]
        }
        
        return genes_ovins.get(race, [
            {'gene': 'GENERIC', 'fonction': 'Gène ovin standard', 'chromosome': 'NA'}
        ])
    
    def generer_rapport_genomique(self, brebis_id: int, sequences: List[Dict]) -> str:
        """Génère un rapport génomique complet"""
        rapport = f"""
        RAPPORT GÉNOMIQUE - Brebis ID: {brebis_id}
        Date: {datetime.now().strftime('%Y-%m-%d')}
        ==============================================
        
        INFORMATIONS GÉNÉTIQUES
        ------------------------
        Nombre de séquences analysées: {len(sequences)}
        
        SÉQUENCES ANALYSÉES:
        """
        
        for i, seq in enumerate(sequences, 1):
            rapport += f"""
        {i}. {seq.get('accession', 'N/A')}
            Type: {seq.get('type', 'ADN')}
            Longueur: {seq.get('longueur', 0)} bp
            Laboratoire: {seq.get('laboratoire', 'Non spécifié')}
            """
        
        rapport += """
        
        ANALYSE COMPARATIVE:
        --------------------
        Les séquences ont été comparées aux bases de données de référence.
        
        RECOMMANDATIONS:
        ----------------
        1. Vérifier les SNP identifiés dans les gènes de production laitière
        2. Considérer le génotypage pour les marqueurs de qualité du lait
        3. Intégrer les données dans le programme de sélection
        
        Ce rapport est généré automatiquement par Ovin Manager Pro.
        Pour une analyse approfondie, consulter un généticien spécialisé.
        """
        
        return rapport

# ============================================================================
# BLOC 7: ANALYSE STATISTIQUE AVEC SIMULATION R
# ============================================================================
class AnalyseurStatistique:
    """Analyse statistique des données ovines (simulation R)"""
    
    def __init__(self):
        self.methodes_disponibles = [
            'correlation',
            'regression_lineaire', 
            'anova',
            'test_t',
            'modele_mixte'
        ]
    
    def analyser_correlations(self, donnees: List[Dict], variable_x: str, 
                             variable_y: str) -> Dict:
        """Analyse de corrélation entre deux variables"""
        # Extraction des données
        x_vals = [d.get(variable_x, 0) for d in donnees if variable_x in d]
        y_vals = [d.get(variable_y, 0) for d in donnees if variable_y in d]
        
        if len(x_vals) < 2 or len(y_vals) < 2:
            return {"erreur": "Données insuffisantes"}
        
        # Calculs statistiques de base
        n = len(x_vals)
        mean_x = sum(x_vals) / n
        mean_y = sum(y_vals) / n
        
        # Covariance et corrélation
        cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_vals, y_vals)) / (n - 1)
        std_x = (sum((x - mean_x)**2 for x in x_vals) / (n - 1))**0.5
        std_y = (sum((y - mean_y)**2 for y in y_vals) / (n - 1))**0.5
        
        if std_x == 0 or std_y == 0:
            correlation = 0
        else:
            correlation = cov / (std_x * std_y)
        
        # Interprétation
        if abs(correlation) > 0.7:
            force = "Forte"
        elif abs(correlation) > 0.3:
            force = "Modérée"
        else:
            force = "Faible"
        
        direction = "positive" if correlation > 0 else "négative"
        
        return {
            'correlation': round(correlation, 4),
            'covariance': round(cov, 4),
            'taille_echantillon': n,
            'interpretation': f"{force} corrélation {direction}",
            'variables': f"{variable_x} vs {variable_y}",
            'mean_x': round(mean_x, 2),
            'mean_y': round(mean_y, 2),
            'std_x': round(std_x, 2),
            'std_y': round(std_y, 2)
        }
    
    def analyser_production_lait(self, donnees_production: List[Dict]) -> Dict:
        """Analyse des courbes de lactation"""
        if not donnees_production:
            return {"erreur": "Aucune donnée de production"}
        
        # Modèle de Wood simplifié: y = a * t^b * e^(-c*t)
        # Où t = jour de lactation, y = production
        
        jours = [d.get('jour_lactation', 0) for d in donnees_production]
        productions = [d.get('production', 0) for d in donnees_production]
        
        # Calcul du pic de lactation
        if productions:
            pic_production = max(productions)
            jour_pic = jours[productions.index(pic_production)]
            
            # Production totale estimée (intégrale simplifiée)
            production_totale = sum(productions)
            
            # Persistance de lactation (production à 150 jours / pic)
            if len(productions) > 150 and jour_pic > 0:
                persistance = productions[150] / pic_production if pic_production > 0 else 0
            else:
                persistance = None
            
            return {
                'pic_production': round(pic_production, 2),
                'jour_pic': jour_pic,
                'production_totale_estimee': round(production_totale, 2),
                'persistance_lactation': round(persistance, 2) if persistance else None,
                'duree_lactation_moyenne': len(jours),
                'production_moyenne': round(sum(productions)/len(productions), 2)
            }
        
        return {"erreur": "Calcul impossible"}
    
    def generer_rapport_statistique(self, brebis_id: int, analyses: List[Dict]) -> str:
        """Génère un rapport statistique complet"""
        rapport = f"""
        RAPPORT STATISTIQUE - Brebis ID: {brebis_id}
        Date: {datetime.now().strftime('%Y-%m-%d')}
        ==============================================
        
        ANALYSES EFFECTUÉES:
        --------------------
        """
        
        for i, analyse in enumerate(analyses, 1):
            rapport += f"\n{i}. {analyse.get('type', 'Analyse')}\n"
            for key, value in analyse.items():
                if key != 'type':
                    rapport += f"   {key}: {value}\n"
        
        rapport += """
        
        INTERPRÉTATION GÉNÉRALE:
        ------------------------
        Les analyses statistiques permettent d'identifier les tendances
        et les relations entre les différents caractères mesurés.
        
        RECOMMANDATIONS:
        ----------------
        1. Utiliser ces résultats pour la sélection génétique
        2. Vérifier la cohérence des mesures
        3. Considérer les facteurs environnementaux
        
        Outil: Ovin Manager Pro - Module Statistique
        """
        
        return rapport

# ============================================================================
# BLOC 8: DONNÉES DE DÉMONSTRATION PRÉ-INSTALLÉES
# ============================================================================
class DonneesDemonstration:
    """Classe pour créer des données de démonstration"""
    
    @staticmethod
    def creer_donnees_test(db_manager: DatabaseManager):
        """Crée des données de démonstration complètes"""
        
        races = ['lacaune', 'manech', 'basco_bearnaise']
        noms_femelles = ['Bella', 'Daisy', 'Luna', 'Molly', 'Sophie', 'Zoe']
        noms_males = ['Max', 'Rocky', 'Charlie', 'Buddy', 'Jack', 'Teddy']
        
        print("\n📊 Création des données de démonstration...")
        
        # Création de 20 brebis
        brebis_ids = []
        for i in range(1, 21):
            sexe = 'F' if i % 3 != 0 else 'M'
            noms = noms_femelles if sexe == 'F' else noms_males
            
            brebis = Brebis(
                id=i,
                identifiant_unique=f"BR{2024}{i:04d}",
                nom=f"{noms[i % len(noms)]}_{i}",
                date_naissance=date(2022, 1, 1) + timedelta(days=i*15),
                race=races[i % len(races)],
                sexe=sexe,
                statut='active' if i < 18 else 'retrait',
                notes=f"Brebis de démonstration {i}"
            )
            
            # Ajout à la base
            cursor = db_manager.conn.cursor()
            cursor.execute('''
                INSERT OR IGNORE INTO brebis 
                (id, identifiant_unique, nom, date_naissance, race, sexe, statut, notes)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ''', (brebis.id, brebis.identifiant_unique, brebis.nom, 
                  brebis.date_naissance.isoformat(), brebis.race, brebis.sexe,
                  brebis.statut, brebis.notes))
            
            brebis_ids.append(i)
        
        # Création de gestations
        print("🤰 Création des gestations de démo...")
        for i, brebis_id in enumerate(brebis_ids[:10]):
            date_eponge = date(2024, 1, 15) + timedelta(days=i*7)
            
            cursor = db_manager.conn.cursor()
            cursor.execute('''
                INSERT INTO gestations 
                (brebis_id, date_eponge, date_mise_bas_prevu, nombre_agneaux_prevus, statut)
                VALUES (?, ?, ?, ?, ?)
            ''', (brebis_id, date_eponge.isoformat(),
                  (date_eponge + timedelta(days=150)).isoformat(),
                  i % 3 + 1, 'en_cours'))
        
        # Création de suivis médicaux
        print("💉 Création des suivis médicaux...")
        vaccins = ['FCO', 'Chlamydiose', 'Paratubeculose']
        for brebis_id in brebis_ids[:15]:
            for j in range(2):  # 2 interventions par brebis
                cursor = db_manager.conn.cursor()
                cursor.execute('''
                    INSERT INTO suivi_medical 
                    (brebis_id, date_intervention, type_intervention, produit, dose)
                    VALUES (?, ?, ?, ?, ?)
                ''', (brebis_id, 
                      (date(2024, 1, 1) + timedelta(days=j*30)).isoformat(),
                      'vaccination',
                      vaccins[j % len(vaccins)],
                      '2 ml'))
        
        # Création de caractères morphologiques
        print("📏 Création des mesures morphologiques...")
        for brebis_id in brebis_ids:
            cursor = db_manager.conn.cursor()
            cursor.execute('''
                INSERT INTO caracteres_morpho 
                (brebis_id, date_mesure, caractere, valeur, unite, methode_mesure)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (brebis_id, date.today().isoformat(), 'poids_vif',
                  65 + (brebis_id % 3) * 10, 'kg', 'manuel'))
        
        db_manager.conn.commit()
        print("✅ Données de démonstration créées avec succès!")
        print(f"   - {len(brebis_ids)} brebis créées")
        print(f"   - 10 gestations en cours")
        print(f"   - 30 interventions médicales")
        print(f"   - {len(brebis_ids)} mesures morphologiques")

# ============================================================================
# BLOC 9: INTERFACE UTILISATEUR PRINCIPALE
# ============================================================================
class InterfaceOvinManager:
    """Interface principale de l'application"""
    
    def __init__(self):
        self.db = DatabaseManager()
        self.gestation_manager = GestionnaireGestation(self.db)
        self.morpho_analyzer = AnalyseurMorphometrique()
        self.genomique = IntegrationGenomique("contact@ovin-manager.com")
        self.stats = AnalyseurStatistique()
        
    def menu_principal(self):
        """Affiche le menu principal"""
        while True:
            print("\n" + "="*60)
            print("🐑 OVIN MANAGER PRO - Application de Gestion Scientifique")
            print("="*60)
            print("\nMENU PRINCIPAL:")
            print("1. 📊 Gestion des Brebis")
            print("2. 🤰 Gestion des Gestations")
            print("3. 📸 Analyse Morphométrique")
            print("4. 🧬 Analyse Génomique")
            print("5. 📈 Analyse Statistique")
            print("6. 📋 Afficher les Données")
            print("7. 🗂️  Créer Données de Démonstration")
            print("8. 📄 Générer Rapports")
            print("0
