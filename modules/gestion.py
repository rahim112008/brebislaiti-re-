"""
Module: Gestion des brebis et gestations
"""

from datetime import datetime, date, timedelta
from typing import Dict, List, Optional
import pandas as pd

class GestionnaireGestation:
    """Gestion des gestations avec prédictions scientifiques"""
    
    DUREES_GESTATION = {
        'lacaune': {'moyenne': 152, 'ecart_type': 2.5},
        'manech': {'moyenne': 150, 'ecart_type': 2.0},
        'basco_bearnaise': {'moyenne': 148, 'ecart_type': 2.2},
        'default': {'moyenne': 150, 'ecart_type': 2.5}
    }
    
    def __init__(self, db_manager):
        self.db = db_manager
    
    def calculer_date_mise_bas(self, date_eponge: date, race: str = 'default') -> date:
        """Calcule la date prévue de mise bas"""
        duree = self.DUREES_GESTATION.get(race, self.DUREES_GESTATION['default'])
        return date_eponge + timedelta(days=int(duree['moyenne']))
    
    def programmer_eponge(self, brebis_id: int, date_eponge: date, race: str) -> Dict:
        """Programme l'introduction d'une éponge et calcule les dates clés"""
        date_mise_bas = self.calculer_date_mise_bas(date_eponge, race)
        date_retrait = date_eponge + timedelta(days=14)
        
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
    
    def get_statistiques_gestation(self) -> Dict:
        """Calcule les statistiques de gestation"""
        query_total = "SELECT COUNT(*) FROM gestations"
        query_en_cours = "SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'"
        query_terminees = "SELECT COUNT(*) FROM gestations WHERE statut = 'termine'"
        query_moyenne = """
            SELECT AVG(nombre_agneaux_nes) 
            FROM gestations 
            WHERE statut = 'termine' AND nombre_agneaux_nes > 0
        """
        
        total = self.db.fetch_one(query_total)[0] or 0
        en_cours = self.db.fetch_one(query_en_cours)[0] or 0
        terminees = self.db.fetch_one(query_terminees)[0] or 0
        moyenne_agneaux = self.db.fetch_one(query_moyenne)[0] or 0
        
        taux_reussite = (terminees / total * 100) if total > 0 else 0
        
        return {
            'total_gestations': total,
            'en_cours': en_cours,
            'terminees': terminees,
            'taux_reussite': taux_reussite,
            'moyenne_agneaux_par_mise_bas': round(moyenne_agneaux, 2)
        }

class DonneesDemonstration:
    """Classe pour créer des données de démonstration"""
    
    @staticmethod
    def creer_donnees_test(db_manager):
        """Crée des données de démonstration complètes"""
        
        races = ['lacaune', 'manech', 'basco_bearnaise']
        noms_femelles = ['Bella', 'Daisy', 'Luna', 'Molly', 'Sophie', 'Zoe']
        noms_males = ['Max', 'Rocky', 'Charlie', 'Buddy', 'Jack', 'Teddy']
        
        print("Création des données de démonstration...")
        
        # Création de 20 brebis
        for i in range(1, 21):
            sexe = 'F' if i % 3 != 0 else 'M'
            noms = noms_femelles if sexe == 'F' else noms_males
            
            identifiant = f"BR{2024}{i:04d}"
            nom = f"{noms[i % len(noms)]}_{i}"
            date_naiss = (date(2022, 1, 1) + timedelta(days=i*15)).isoformat()
            race = races[i % len(races)]
            statut = 'active' if i < 18 else 'retrait'
            
            query = '''
                INSERT OR IGNORE INTO brebis 
                (identifiant_unique, nom, date_naissance, race, sexe, statut, notes)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            '''
            
            db_manager.execute_query(query, (
                identifiant, nom, date_naiss, race, sexe, statut,
                f"Brebis de démonstration {i}"
            ))
        
        # Création de gestations
        for i in range(1, 11):
            date_eponge = (date(2024, 1, 15) + timedelta(days=i*7)).isoformat()
            date_mise_bas = (date(2024, 6, 15) + timedelta(days=i*7)).isoformat()
            
            query = '''
                INSERT INTO gestations 
                (brebis_id, date_eponge, date_mise_bas_prevu, nombre_agneaux_prevus, statut)
                VALUES (?, ?, ?, ?, ?)
            '''
            
            db_manager.execute_query(query, (i, date_eponge, date_mise_bas, i % 3 + 1, 'en_cours'))
        
        print("✅ Données de démonstration créées")
