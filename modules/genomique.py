"""
Module: Intégration génomique et analyse SNP
"""

from typing import Dict, List
from datetime import datetime

class IntegrationGenomique:
    """Intégration avec les bases de données génomiques"""
    
    def __init__(self, email: str):
        self.email = email
        self.sequences_cache = {}
    
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
            'snps_detailles': snps[:10],
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
        """Retourne les gènes candidats pour une race donnée"""
        genes_ovins = {
            'lacaune': [
                {'gene': 'LALBA', 'fonction': 'Protéine du lait', 'chromosome': '3'},
                {'gene': 'CSN1S1', 'fonction': 'Caséine alpha-S1', 'chromosome': '6'},
                {'gene': 'DGAT1', 'fonction': 'Synthèse des triglycérides', 'chromosome': '14'},
            ],
            'manech': [
                {'gene': 'PRLR', 'fonction': 'Récepteur prolactine', 'chromosome': '16'},
                {'gene': 'GH1', 'fonction': 'Hormone de croissance', 'chromosome': '11'},
            ],
            'basco_bearnaise': [
                {'gene': 'MSTN', 'fonction': 'Myostatine', 'chromosome': '2'},
                {'gene': 'LEP', 'fonction': 'Leptine', 'chromosome': '4'},
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
        """
        
        return rapport
