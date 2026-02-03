"""
Module: Analyse statistique
"""

from typing import Dict, List
import numpy as np

class AnalyseurStatistique:
    """Analyse statistique des données ovines"""
    
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
        
        # Calculs statistiques
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
        
        jours = [d.get('jour_lactation', 0) for d in donnees_production]
        productions = [d.get('production', 0) for d in donnees_production]
        
        if productions:
            pic_production = max(productions)
            jour_pic = jours[productions.index(pic_production)]
            
            production_totale = sum(productions)
            
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
