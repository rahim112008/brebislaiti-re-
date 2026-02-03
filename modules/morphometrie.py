"""
Module: Analyse morphométrique par smartphone
"""

from typing import Dict
from datetime import datetime

class AnalyseurMorphometrique:
    """Analyse morphométrique à partir de photos smartphone"""
    
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
            'tour_poitrine': 'cm',
            'poids_estime': 'kg'
        }
    
    def analyser_photo(self, image_path: str, objet_reference_pixels: int, 
                      taille_reelle_objet: float, race: str = None) -> Dict:
        """
        Analyse une photo pour extraire des mesures morphométriques
        """
        # Facteur de conversion pixels -> cm
        facteur_conversion = taille_reelle_objet / objet_reference_pixels
        
        # Simulation d'analyse d'image
        mesures = {}
        
        # Longueur du corps
        longueur_pixels = objet_reference_pixels * 8
        mesures['longueur_corps'] = round(longueur_pixels * facteur_conversion, 2)
        
        # Hauteur au garrot
        mesures['hauteur_garrot'] = round(mesures['longueur_corps'] * 0.85, 2)
        
        # Tour de poitrine
        mesures['tour_poitrine'] = round(mesures['longueur_corps'] * 1.2, 2)
        
        # Poids estimé
        poids_estime = (mesures['tour_poitrine']**2 * mesures['longueur_corps']) / 30000
        mesures['poids_estime'] = round(poids_estime, 2)
        
        # Analyse comparative si race spécifiée
        if race and race in self.REFERENCES_RACES:
            reference = self.REFERENCES_RACES[race]
            analyses = []
            
            if 'poids_adulte_femelle' in reference:
                min_ref, max_ref = reference['poids_adulte_femelle']
                if min_ref <= mesures['poids_estime'] <= max_ref:
                    analyses.append(f"Poids dans la norme pour {race}")
                else:
                    analyses.append(f"Poids hors norme pour {race} ({min_ref}-{max_ref} kg)")
            
            return {
                'mesures': mesures,
                'unites': self.mesures_standard,
                'facteur_conversion': facteur_conversion,
                'analyse_comparative': analyses,
                'date_analyse': datetime.now().isoformat(),
                'precision': 'estimative'
            }
        
        return {
            'mesures': mesures,
            'unites': self.mesures_standard,
            'facteur_conversion': facteur_conversion,
            'date_analyse': datetime.now().isoformat()
        }
