"""
OVIN MANAGER PRO - Version Phénotypique Complète avec Module Génétique Avancé
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
from typing import Dict, List, Tuple, Optional
import plotly.figure_factory as ff
import requests
import io
import base64
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import networkx as nx

# ========== CONFIGURATION ==========

st.set_page_config(
    page_title="Ovin Manager Pro - Génétique Avancée",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personnalisé
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #2E7D32;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.8rem;
        color: #388E3C;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 10px;
        padding: 15px;
        border-left: 5px solid #28a745;
    }
    .module-card {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        border-radius: 10px;
        padding: 20px;
        margin: 10px 0;
    }
    .dna-sequence {
        font-family: 'Courier New', monospace;
        background-color: #f0f0f0;
        padding: 10px;
        border-radius: 5px;
        margin: 5px 0;
        font-size: 0.9em;
        letter-spacing: 1px;
    }
    .gene-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

# ========== INITIALISATION BASE DE DONNÉES ==========

def init_database():
    """Initialise la base de données SQLite complète"""
    conn = sqlite3.connect('ovin_manager_genetic.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table des brebis
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT,
            date_naissance DATE,
            race TEXT,
            sexe TEXT,
            statut TEXT DEFAULT 'active',
            poids FLOAT,
            notes TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Table gestations
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_eponge DATE,
            date_mise_bas_prevu DATE,
            nombre_agneaux_prevus INTEGER DEFAULT 1,
            statut TEXT DEFAULT 'en_cours',
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
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
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
    ''')
    
    # Table données génomiques avancées
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS donnees_genomiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            gene_nom TEXT,
            sequence_adn TEXT,
            chromosome TEXT,
            position_start INTEGER,
            position_end INTEGER,
            type_mutation TEXT,
            allele1 TEXT,
            allele2 TEXT,
            genotype TEXT,
            frequence_allele FLOAT,
            effet_phenotype TEXT,
            qualite_score INTEGER,
            date_analyse DATE,
            source_db TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
    ''')
    
    # Table marqueurs SNP
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS snp_marqueurs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            rs_id TEXT UNIQUE,
            chromosome TEXT,
            position INTEGER,
            allele_reference TEXT,
            allele_alternatif TEXT,
            gene_associe TEXT,
            fonction TEXT,
            impact TEXT,
            frequence_maf FLOAT,
            heritabilite FLOAT,
            qtl_associe TEXT,
            date_ajout DATE
        )
    ''')
    
    # Table QTL (Quantitative Trait Loci)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS qtl_ovins (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nom_qtl TEXT,
            chromosome TEXT,
            position_start INTEGER,
            position_end INTEGER,
            caractere_etudie TEXT,
            lods_score FLOAT,
            variance_expliquee FLOAT,
            race_etudiee TEXT,
            publication_reference TEXT,
            genes_candidats TEXT
        )
    ''')
    
    # Table analyses génétiques
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS analyses_genetiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            type_analyse TEXT,
            resultats_json TEXT,
            score_genetique FLOAT,
            recommendations TEXT,
            date_analyse DATE,
            analyse_par TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
    ''')
    
    # Table pedigrees
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pedigrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            animal_id TEXT UNIQUE,
            pere_id TEXT,
            mere_id TEXT,
            race TEXT,
            generation INTEGER,
            coefficient_consanguinite FLOAT,
            date_naissance DATE,
            FOREIGN KEY (animal_id) REFERENCES brebis(identifiant_unique)
        )
    ''')
    
    conn.commit()
    return conn

# Connexion à la base de données
conn = init_database()

# ========== MODULE GÉNÉTIQUE AVANCÉ ==========

class GeneticAnalyzer:
    """Analyseur génétique avancé pour professionnels"""
    
    @staticmethod
    def sequence_analyzer(sequence: str) -> Dict:
        """Analyse approfondie d'une séquence ADN"""
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        
        # Composition nucléotidique
        composition = {
            'A': seq.count('A'),
            'T': seq.count('T'),
            'C': seq.count('C'),
            'G': seq.count('G'),
            'N': seq.count('N') + seq.count('X')
        }
        
        total = sum(composition.values())
        
        # Calculs avancés
        gc_content = ((composition['G'] + composition['C']) / total * 100) if total > 0 else 0
        at_content = 100 - gc_content
        gc_skew = (composition['G'] - composition['C']) / (composition['G'] + composition['C']) if (composition['G'] + composition['C']) > 0 else 0
        
        # Recherche de motifs
        motifs = {
            'start_codon': seq.count('ATG'),
            'stop_codons': seq.count('TAA') + seq.count('TAG') + seq.count('TGA'),
            'cpgi': GeneticAnalyzer._find_cpg_islands(seq),
            'repeats': GeneticAnalyzer._find_repeats(seq),
            'restriction_sites': GeneticAnalyzer._find_restriction_sites(seq)
        }
        
        # Prédiction de caractéristiques
        prediction = {
            'is_coding': GeneticAnalyzer._predict_coding_potential(seq),
            'melting_temp': GeneticAnalyzer._calculate_tm(seq),
            'molecular_weight': GeneticAnalyzer._calculate_mw(seq),
            'secondary_structure': GeneticAnalyzer._predict_secondary_structure(seq)
        }
        
        return {
            'longueur': total,
            'composition': composition,
            'pourcentages': {
                'GC': round(gc_content, 2),
                'AT': round(at_content, 2),
                'GC_skew': round(gc_skew, 3)
            },
            'motifs': motifs,
            'prediction': prediction,
            'checksum': hash(seq) % 10000
        }
    
    @staticmethod
    def _find_cpg_islands(seq: str, window=200, threshold=0.6) -> List[Dict]:
        """Détecte les îlots CpG"""
        islands = []
        for i in range(0, len(seq) - window + 1, window//2):
            window_seq = seq[i:i+window]
            cpg_count = window_seq.count('CG')
            gc_content = (window_seq.count('G') + window_seq.count('C')) / len(window_seq)
            
            if cpg_count > 0 and gc_content > threshold:
                islands.append({
                    'start': i,
                    'end': i + window,
                    'cpg_count': cpg_count,
                    'gc_content': round(gc_content, 3)
                })
        return islands[:5]  # Retourne les 5 premiers
    
    @staticmethod
    def _find_repeats(seq: str) -> Dict:
        """Trouve les séquences répétées"""
        repeats = {}
        for length in [2, 3, 4]:
            repeat_dict = {}
            for i in range(0, len(seq) - length + 1):
                motif = seq[i:i+length]
                if motif in repeat_dict:
                    repeat_dict[motif] += 1
                else:
                    repeat_dict[motif] = 1
            
            # Garder les motifs les plus fréquents
            top_motifs = sorted(repeat_dict.items(), key=lambda x: x[1], reverse=True)[:3]
            repeats[f'{length}mer'] = top_motifs
        
        return repeats
    
    @staticmethod
    def _find_restriction_sites(seq: str) -> List[Dict]:
        """Trouve les sites de restriction courants"""
        enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XbaI': 'TCTAGA',
            'NotI': 'GCGGCCGC'
        }
        
        sites = []
        for enzyme, site in enzymes.items():
            positions = [i for i in range(len(seq) - len(site) + 1) if seq[i:i+len(site)] == site]
            if positions:
                sites.append({
                    'enzyme': enzyme,
                    'site': site,
                    'positions': positions,
                    'count': len(positions)
                })
        
        return sites
    
    @staticmethod
    def _predict_coding_potential(seq: str) -> Dict:
        """Prédit le potentiel de codage"""
        # Algorithme simplifié basé sur la périodicité
        frames = []
        for frame in range(3):
            codons = [seq[i:i+3] for i in range(frame, len(seq)-2, 3)]
            stop_count = sum(1 for codon in codons if codon in ['TAA', 'TAG', 'TGA'])
            frames.append({
                'frame': frame + 1,
                'stop_codons': stop_count,
                'coding_score': max(0, 1 - (stop_count / max(1, len(codons)/10)))
            })
        
        best_frame = max(frames, key=lambda x: x['coding_score'])
        return {
            'frames': frames,
            'best_frame': best_frame,
            'is_likely_coding': best_frame['coding_score'] > 0.7
        }
    
    @staticmethod
    def _calculate_tm(seq: str) -> float:
        """Calcule la température de fusion (formule Wallace)"""
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        return 2 * at_count + 4 * gc_count
    
    @staticmethod
    def _calculate_mw(seq: str) -> float:
        """Calcule le poids moléculaire"""
        weights = {'A': 313.21, 'T': 304.2, 'C': 289.18, 'G': 329.21}
        total = sum(weights.get(base, 300) for base in seq)
        return total / 1000  # en kDa
    
    @staticmethod
    def _predict_secondary_structure(seq: str) -> Dict:
        """Prédit la structure secondaire simplifiée"""
        # Prédiction basée sur la composition
        gc_content = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
        
        if gc_content > 0.6:
            structure = "Fortement structuré (GC-rich)"
        elif gc_content > 0.4:
            structure = "Modérément structuré"
        else:
            structure = "Peu structuré (AT-rich)"
        
        return {
            'predicted_structure': structure,
            'gc_content': gc_content,
            'stem_loop_potential': round(gc_content * 100, 1)
        }
    
    @staticmethod
    def alignment_analyzer(seq1: str, seq2: str) -> Dict:
        """Alignement de séquences avec analyse détaillée"""
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        alignments = aligner.align(seq1, seq2)
        
        if alignments:
            alignment = alignments[0]
            alignment_str = str(alignment).split('\n')
            
            # Statistiques d'alignement
            matches = sum(1 for a, b in zip(alignment_str[0], alignment_str[2]) if a == b)
            gaps = alignment_str[1].count('-')
            total = len(alignment_str[0])
            
            return {
                'score': alignment.score,
                'identity': round(matches / total * 100, 2),
                'gaps': gaps,
                'gap_percentage': round(gaps / total * 100, 2),
                'length': total,
                'alignment': alignment_str[:3],
                'coverage': round(len(seq1) / total * 100, 2)
            }
        
        return {'error': 'Alignement impossible'}
    
    @staticmethod
    def pedigree_analyzer(pedigree_data: List[Tuple]) -> Dict:
        """Analyse de pedigree avec calculs de consanguinité"""
        # Construction du graphe de pedigree
        G = nx.DiGraph()
        
        for animal, sire, dam in pedigree_data:
            G.add_node(animal)
            if sire:
                G.add_edge(sire, animal)
            if dam:
                G.add_edge(dam, animal)
        
        # Calculs avancés
        inbreeding_coeffs = {}
        for node in G.nodes():
            # Coefficient de consanguinité simplifié
            ancestors = list(nx.ancestors(G, node))
            if len(ancestors) > 1:
                # Calcul basique
                inbreeding_coeffs[node] = round(1 / (2 ** len(ancestors)), 4)
            else:
                inbreeding_coeffs[node] = 0.0
        
        return {
            'total_animals': len(G.nodes()),
            'total_relations': len(G.edges()),
            'inbreeding_coefficients': inbreeding_coeffs,
            'average_inbreeding': round(np.mean(list(inbreeding_coeffs.values())), 4),
            'generations': GeneticAnalyzer._count_generations(G),
            'founder_animals': GeneticAnalyzer._find_founders(G)
        }
    
    @staticmethod
    def _count_generations(G):
        """Compte les générations dans le pedigree"""
        generations = {}
        for node in G.nodes():
            depth = len(list(nx.shortest_path_length(G, node).values()))
            generations[node] = depth
        return generations
    
    @staticmethod
    def _find_founders(G):
        """Trouve les animaux fondateurs"""
        return [node for node in G.nodes() if G.in_degree(node) == 0]

class NCBIIntegration:
    """Intégration avec les bases de données NCBI"""
    
    @staticmethod
    def search_ncbi(query: str, db: str = "nuccore", retmax: int = 10) -> List[Dict]:
        """Recherche dans NCBI (version simulée pour l'exemple)"""
        
        # Données simulées pour démonstration
        mock_results = [
            {
                'id': 'NC_019458.2',
                'title': 'Ovis aries breed Romanov chromosome 1, whole genome shotgun sequence',
                'species': 'Ovis aries',
                'length': 275612895,
                'date': '2013/12/20',
                'features': ['genes', 'CDS', 'mRNA', 'tRNA']
            },
            {
                'id': 'NC_019459.2',
                'title': 'Ovis aries breed Romanov chromosome 2',
                'species': 'Ovis aries',
                'length': 248993846,
                'date': '2013/12/20',
                'features': ['genes', 'repeats', 'SNPs']
            },
            {
                'id': 'XM_004005000.3',
                'title': 'Ovis aries growth differentiation factor 8 (GDF8), mRNA',
                'species': 'Ovis aries',
                'length': 2856,
                'date': '2022/05/15',
                'features': ['CDS', 'exons', 'UTR']
            },
            {
                'id': 'NM_001009394.1',
                'title': 'Ovis aries myostatin (MSTN), mRNA',
                'species': 'Ovis aries',
                'length': 1128,
                'date': '2006/04/02',
                'features': ['coding', 'polypeptide']
            }
        ]
        
        # Filtrer par requête
        filtered = [res for res in mock_results if query.lower() in str(res).lower()]
        return filtered[:retmax]
    
    @staticmethod
    def get_gene_info(gene_id: str) -> Dict:
        """Récupère les informations d'un gène"""
        gene_database = {
            'MSTN': {
                'nom': 'Myostatine',
                'synonymes': ['GDF8', 'Growth Differentiation Factor 8'],
                'chromosome': '2',
                'position': '6254871-6265123',
                'fonction': 'Régulateur négatif de la croissance musculaire',
                'phenotypes': ['Hypertrophie musculaire', 'Double-muscling'],
                'mutations_connues': ['g.6723G>A', 'c.939G>A'],
                'heritabilite': 0.85
            },
            'PRNP': {
                'nom': 'Prion Protein',
                'synonymes': ['PrP'],
                'chromosome': '13',
                'position': '42316543-42328976',
                'fonction': 'Protéine prion, susceptibilité aux encéphalopathies',
                'phenotypes': ['Résistance/tolérance à la tremblante'],
                'mutations_connues': ['codon 136', 'codon 154', 'codon 171'],
                'heritabilite': 0.92
            },
            'DGAT1': {
                'nom': 'Diacylglycerol O-Acyltransferase 1',
                'chromosome': '14',
                'position': '21894765-21912345',
                'fonction': 'Synthèse des triglycérides',
                'phenotypes': ['Teneur en matière grasse du lait'],
                'mutations_connues': ['K232A'],
                'heritabilite': 0.45
            }
        }
        
        return gene_database.get(gene_id, {'error': 'Gène non trouvé dans la base'})
    
    @staticmethod
    def fetch_sequence(accession: str) -> Dict:
        """Récupère une séquence (simulée)"""
        sequences = {
            'NC_019458.2': 'ATCG' * 1000,
            'XM_004005000.3': 'ATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGCCTGGGCCGCATGG' * 50,
            'NM_001009394.1': 'ATG' + 'GCT' * 375  # Séquence MSTN simplifiée
        }
        
        seq = sequences.get(accession, '')
        return {
            'accession': accession,
            'sequence': seq,
            'length': len(seq),
            'gc_content': round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2) if seq else 0
        }

class PopulationGenetics:
    """Analyses de génétique des populations"""
    
    @staticmethod
    def hardy_weinberg(genotypes: List[str]) -> Dict:
        """Test d'équilibre de Hardy-Weinberg"""
        from collections import Counter
        
        counts = Counter(genotypes)
        total = sum(counts.values())
        
        # Calcul des fréquences alléliques
        alleles = []
        for genotype in genotypes:
            alleles.extend(list(genotype))
        
        allele_counts = Counter(alleles)
        allele_freq = {allele: count/(total*2) for allele, count in allele_counts.items()}
        
        # Fréquences attendues
        expected = {}
        for a1 in allele_freq:
            for a2 in allele_freq:
                genotype = ''.join(sorted([a1, a2]))
                freq = allele_freq[a1] * allele_freq[a2]
                if a1 != a2:
                    freq *= 2
                expected[genotype] = freq * total
        
        # Test du chi²
        chi2 = 0
        for genotype in set(list(counts.keys()) + list(expected.keys())):
            obs = counts.get(genotype, 0)
            exp = expected.get(genotype, 0)
            if exp > 0:
                chi2 += ((obs - exp) ** 2) / exp
        
        return {
            'allele_frequencies': allele_freq,
            'observed': counts,
            'expected': {k: round(v, 2) for k, v in expected.items()},
            'chi_squared': round(chi2, 4),
            'p_value': round(stats.chi2.sf(chi2, df=1), 4),
            'in_hardy_weinberg': stats.chi2.sf(chi2, df=1) > 0.05
        }
    
    @staticmethod
    def genetic_diversity(genotypes: List[str]) -> Dict:
        """Mesures de diversité génétique"""
        # Nombre d'allèles
        alleles = set()
        for genotype in genotypes:
            alleles.update(list(genotype))
        
        # Hétérozygotie observée et attendue
        het_obs = sum(1 for g in genotypes if len(set(g)) > 1) / len(genotypes)
        
        # Fréquences alléliques
        allele_counts = {}
        total_alleles = len(genotypes) * 2
        for genotype in genotypes:
            for allele in genotype:
                allele_counts[allele] = allele_counts.get(allele, 0) + 1
        
        allele_freq = {a: c/total_alleles for a, c in allele_counts.items()}
        
        # Hétérozygotie attendue
        het_exp = 1 - sum(f**2 for f in allele_freq.values())
        
        # F-statistiques
        fis = 1 - (het_obs / het_exp) if het_exp > 0 else 0
        
        return {
            'allele_count': len(alleles),
            'heterozygosity_observed': round(het_obs, 4),
            'heterozygosity_expected': round(het_exp, 4),
            'fis_inbreeding': round(fis, 4),
            'allele_frequencies': allele_freq,
            'shannon_index': round(-sum(f * np.log(f) for f in allele_freq.values()), 4)
        }
    
    @staticmethod
    def pca_analysis(genotype_matrix: np.ndarray) -> Dict:
        """Analyse en composantes principales"""
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(genotype_matrix)
        
        pca = PCA(n_components=3)
        components = pca.fit_transform(scaled_data)
        
        return {
            'explained_variance': pca.explained_variance_ratio_.tolist(),
            'components': components.tolist(),
            'total_variance': sum(pca.explained_variance_ratio_),
            'loadings': pca.components_.tolist()
        }

# ========== PAGES DE L'APPLICATION ==========

def afficher_genetique_avancee():
    """Page principale de génétique avancée"""
    
    st.markdown('<h1 class="main-header">🧬 Module de Génétique Avancée</h1>', unsafe_allow_html=True)
    st.markdown("*Pour généticiens et chercheurs - Analyses professionnelles*")
    
    # Menu de navigation génétique
    genetic_tabs = st.tabs([
        "🧬 Analyse Séquences", 
        "🔍 Recherche NCBI", 
        "📊 Génétique Pop.", 
        "🧮 SNP & QTL",
        "🌳 Pedigrees",
        "📈 GWAS",
        "💾 Import/Export"
    ])
    
    # Tab 1: Analyse de séquences
    with genetic_tabs[0]:
        st.markdown("### Analyse Avancée de Séquences ADN")
        
        col_seq1, col_seq2 = st.columns([2, 1])
        
        with col_seq1:
            sequence_input = st.text_area(
                "Collez votre séquence ADN (FASTA ou format brut):",
                height=250,
                placeholder=">Sequence_Ovis_aries_gene_X\nATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGC..."
            )
            
            # Options d'analyse
            analysis_options = st.multiselect(
                "Types d'analyse:",
                ["Composition", "Motifs", "Structure secondaire", "Sites de restriction", 
                 "Potentiel codant", "Alignement", "Traduction"]
            )
        
        with col_seq2:
            st.markdown("#### Bases de données de référence")
            
            reference_db = st.selectbox(
                "Séquence de référence:",
                ["Ovis aries (GCF_000298735.2)", "Bos taurus (GCF_002263795.2)", 
                 "Homo sapiens (GCF_000001405.40)", "Aucune"]
            )
            
            st.markdown("---")
            st.markdown("#### Outils")
            
            if st.button("🧪 Analyser la séquence", type="primary"):
                if sequence_input:
                    with st.spinner("Analyse en cours..."):
                        # Nettoyer la séquence
                        lines = sequence_input.strip().split('\n')
                        sequence = ''.join([line for line in lines if not line.startswith('>')])
                        sequence = sequence.upper().replace(" ", "").replace("\n", "")
                        
                        if len(sequence) > 0:
                            # Analyse complète
                            results = GeneticAnalyzer.sequence_analyzer(sequence)
                            
                            # Afficher les résultats
                            st.success(f"✅ Analyse terminée! Séquence de {results['longueur']} bp")
                            
                            # Métriques principales
                            col_m1, col_m2, col_m3, col_m4 = st.columns(4)
                            
                            with col_m1:
                                st.metric("Longueur", f"{results['longueur']} bp")
                            with col_m2:
                                st.metric("% GC", f"{results['pourcentages']['GC']}%")
                            with col_m3:
                                st.metric("Poids moléculaire", f"{results['prediction']['molecular_weight']:.1f} kDa")
                            with col_m4:
                                st.metric("Température de fusion", f"{results['prediction']['melting_temp']:.1f}°C")
                            
                            # Onglets détaillés
                            result_tabs = st.tabs(["Composition", "Motifs", "Structure", "Rapport"])
                            
                            with result_tabs[0]:
                                # Graphique de composition
                                fig = px.pie(
                                    values=list(results['composition'].values()),
                                    names=list(results['composition'].keys()),
                                    title="Composition nucléotidique"
                                )
                                st.plotly_chart(fig, use_container_width=True)
                                
                                st.write("**Détails:**")
                                st.json(results['composition'])
                            
                            with result_tabs[1]:
                                if results['motifs']:
                                    st.write("**Sites de restriction trouvés:**")
                                    for site in results['motifs']['restriction_sites']:
                                        st.info(f"**{site['enzyme']}**: {site['site']} à {len(site['positions'])} position(s)")
                                    
                                    st.write("**Répétitions:**")
                                    for repeat_type, motifs in results['motifs']['repeats'].items():
                                        if motifs:
                                            st.write(f"**{repeat_type}**:")
                                            for motif, count in motifs:
                                                st.write(f"  - {motif}: {count} occurrences")
                                else:
                                    st.info("Aucun motif spécifique détecté")
                            
                            with result_tabs[2]:
                                st.write("**Prédiction de structure secondaire:**")
                                st.success(results['prediction']['secondary_structure']['predicted_structure'])
                                
                                # Graphique GC skew
                                if len(sequence) > 100:
                                    window_size = min(100, len(sequence)//10)
                                    gc_skews = []
                                    positions = []
                                    
                                    for i in range(0, len(sequence)-window_size+1, window_size//2):
                                        window = sequence[i:i+window_size]
                                        gc = window.count('G') + window.count('C')
                                        at = window.count('A') + window.count('T')
                                        skew = (window.count('G') - window.count('C')) / max(1, gc)
                                        gc_skews.append(skew)
                                        positions.append(i)
                                    
                                    fig = go.Figure()
                                    fig.add_trace(go.Scatter(x=positions, y=gc_skews, 
                                                           mode='lines', name='GC Skew'))
                                    fig.update_layout(title="GC Skew le long de la séquence",
                                                    xaxis_title="Position",
                                                    yaxis_title="GC Skew")
                                    st.plotly_chart(fig, use_container_width=True)
                            
                            with result_tabs[3]:
                                # Rapport téléchargeable
                                report = {
                                    "analyse_date": datetime.now().isoformat(),
                                    "sequence_length": results['longueur'],
                                    "gc_content": results['pourcentages']['GC'],
                                    "analysis_results": results
                                }
                                
                                st.download_button(
                                    label="📥 Télécharger le rapport complet (JSON)",
                                    data=json.dumps(report, indent=2),
                                    file_name=f"genetic_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                                    mime="application/json"
                                )
                                
                                st.code(json.dumps(results, indent=2), language='json')
                        else:
                            st.error("Séquence vide ou invalide")
                else:
                    st.warning("Veuillez entrer une séquence ADN")
            
            # Outil d'alignement
            st.markdown("---")
            st.markdown("#### Alignement de séquences")
            seq1 = st.text_input("Séquence 1:", placeholder="ATCG...")
            seq2 = st.text_input("Séquence 2:", placeholder="ATCG...")
            
            if st.button("⚡ Aligner"):
                if seq1 and seq2:
                    alignment = GeneticAnalyzer.alignment_analyzer(seq1, seq2)
                    if 'error' not in alignment:
                        st.metric("Identité", f"{alignment['identity']}%")
                        st.metric("Score", f"{alignment['score']}")
                        
                        # Afficher l'alignement
                        st.write("**Alignement:**")
                        for line in alignment['alignment']:
                            st.code(line)
    
    # Tab 2: Recherche NCBI
    with genetic_tabs[1]:
        st.markdown("### 🔍 Recherche dans les Bases de Données NCBI")
        
        col_ncbi1, col_ncbi2 = st.columns([3, 1])
        
        with col_ncbi1:
            query = st.text_input("Terme de recherche:", 
                                 value="Ovis aries MSTN", 
                                 placeholder="Ex: Ovis aries myostatin gene")
            
            db_options = st.multiselect(
                "Bases de données:",
                ["nuccore", "nucleotide", "gene", "protein", "genome", "snp"],
                default=["nuccore", "gene"]
            )
        
        with col_ncbi2:
            max_results = st.slider("Nombre de résultats", 5, 50, 20)
            sort_by = st.selectbox("Trier par:", ["Pertinence", "Date", "Longueur"])
        
        if st.button("🔬 Rechercher NCBI", type="primary"):
            with st.spinner("Recherche en cours..."):
                results = NCBIIntegration.search_ncbi(query, retmax=max_results)
                
                if results:
                    st.success(f"✅ {len(results)} résultat(s) trouvé(s)")
                    
                    for i, result in enumerate(results, 1):
                        with st.expander(f"Résultat {i}: {result['title'][:100]}..."):
                            col_res1, col_res2 = st.columns([2, 1])
                            
                            with col_res1:
                                st.write(f"**ID:** `{result['id']}`")
                                st.write(f"**Espèce:** {result['species']}")
                                st.write(f"**Longueur:** {result['length']:,} bp")
                                st.write(f"**Date:** {result['date']}")
                                st.write(f"**Caractéristiques:** {', '.join(result['features'])}")
                            
                            with col_res2:
                                # Boutons d'action
                                if st.button(f"📥 Récupérer séquence", key=f"fetch_{i}"):
                                    seq_data = NCBIIntegration.fetch_sequence(result['id'])
                                    st.session_state[f'sequence_{i}'] = seq_data
                                    st.success(f"Séquence de {seq_data['length']} bp récupérée!")
                                
                                if st.button(f"💾 Sauvegarder localement", key=f"save_{i}"):
                                    # Sauvegarde dans la base
                                    cursor = conn.cursor()
                                    try:
                                        cursor.execute('''
                                            INSERT INTO donnees_genomiques 
                                            (gene_nom, sequence_adn, chromosome, date_analyse, source_db)
                                            VALUES (?, ?, ?, ?, ?)
                                        ''', (
                                            result['title'].split()[0],
                                            NCBIIntegration.fetch_sequence(result['id'])['sequence'][:1000],
                                            'Unknown',
                                            date.today().isoformat(),
                                            'NCBI'
                                        ))
                                        conn.commit()
                                        st.success("Séquence sauvegardée!")
                                    except Exception as e:
                                        st.error(f"Erreur: {e}")
                
                # Recherche de gènes spécifiques
                st.markdown("### 🧬 Recherche de Gènes Ovin")
                
                gene_search = st.text_input("Nom du gène:", placeholder="Ex: MSTN, PRNP, DGAT1...")
                
                if gene_search:
                    gene_info = NCBIIntegration.get_gene_info(gene_search.upper())
                    
                    if 'error' not in gene_info:
                        st.markdown(f"#### **{gene_search.upper()}** - {gene_info['nom']}")
                        
                        col_gene1, col_gene2 = st.columns(2)
                        
                        with col_gene1:
                            st.write(f"**Chromosome:** {gene_info['chromosome']}")
                            st.write(f"**Position:** {gene_info['position']}")
                            st.write(f"**Fonction:** {gene_info['fonction']}")
                            st.write(f"**Héritabilité:** {gene_info['heritabilite']}")
                        
                        with col_gene2:
                            st.write("**Phénotypes associés:**")
                            for pheno in gene_info['phenotypes']:
                                st.write(f"- {pheno}")
                            
                            st.write("**Mutations connues:**")
                            for mut in gene_info['mutations_connues']:
                                st.write(f"- {mut}")
    
    # Tab 3: Génétique des populations
    with genetic_tabs[2]:
        st.markdown("### 📊 Génétique des Populations")
        
        pop_tabs = st.tabs(["Hardy-Weinberg", "Diversité", "PCA", "Structure"])
        
        with pop_tabs[0]:
            st.markdown("#### Test d'Équilibre Hardy-Weinberg")
            
            # Entrée de génotypes
            genotypes_input = st.text_area(
                "Entrez les génotypes (un par ligne, ex: AA, AB, BB):",
                height=150,
                placeholder="AA\nAB\nBB\nAA\nAB\nAA\nBB\nAB\nAA"
            )
            
            if genotypes_input:
                genotypes = [g.strip().upper() for g in genotypes_input.split('\n') if g.strip()]
                
                if st.button("📊 Calculer H-W"):
                    hw_results = PopulationGenetics.hardy_weinberg(genotypes)
                    
                    col_hw1, col_hw2, col_hw3 = st.columns(3)
                    
                    with col_hw1:
                        st.metric("χ²", f"{hw_results['chi_squared']:.4f}")
                    with col_hw2:
                        st.metric("p-value", f"{hw_results['p_value']:.4f}")
                    with col_hw3:
                        status = "✅ Équilibre" if hw_results['in_hardy_weinberg'] else "⚠️ Déséquilibre"
                        st.metric("Statut", status)
                    
                    # Graphique observé vs attendu
                    genotypes_list = list(set(genotypes))
                    observed = [hw_results['observed'].get(g, 0) for g in genotypes_list]
                    expected = [hw_results['expected'].get(g, 0) for g in genotypes_list]
                    
                    fig = go.Figure(data=[
                        go.Bar(name='Observé', x=genotypes_list, y=observed),
                        go.Bar(name='Attendu', x=genotypes_list, y=expected)
                    ])
                    
                    fig.update_layout(
                        title="Distribution des génotypes - Observé vs Attendu",
                        barmode='group'
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
        
        with pop_tabs[1]:
            st.markdown("#### Analyse de Diversité Génétique")
            
            if 'genotypes' in locals() and genotypes:
                diversity = PopulationGenetics.genetic_diversity(genotypes)
                
                col_div1, col_div2, col_div3 = st.columns(3)
                
                with col_div1:
                    st.metric("Nombre d'allèles", diversity['allele_count'])
                with col_div2:
                    st.metric("Hétérozygotie observée", f"{diversity['heterozygosity_observed']:.4f}")
                with col_div3:
                    st.metric("Fis", f"{diversity['fis_inbreeding']:.4f}")
                
                # Graphique des fréquences alléliques
                fig = px.bar(
                    x=list(diversity['allele_frequencies'].keys()),
                    y=list(diversity['allele_frequencies'].values()),
                    title="Fréquences alléliques"
                )
                st.plotly_chart(fig, use_container_width=True)
        
        with pop_tabs[2]:
            st.markdown("#### Analyse en Composantes Principales (PCA)")
            
            # Données d'exemple
            st.info("Chargement de données d'exemple...")
            
            # Génération de données simulées
            np.random.seed(42)
            n_animals = 50
            n_snps = 100
            
            genotype_matrix = np.random.choice([0, 1, 2], size=(n_animals, n_snps))
            races = np.random.choice(['Ouled_Djellal', 'Razè', 'Hamra', 'Dman'], size=n_animals)
            
            if st.button("📈 Lancer PCA"):
                pca_results = PopulationGenetics.pca_analysis(genotype_matrix)
                
                # Graphique PCA
                df_pca = pd.DataFrame({
                    'PC1': [x[0] for x in pca_results['components']],
                    'PC2': [x[1] for x in pca_results['components']],
                    'Race': races
                })
                
                fig = px.scatter(df_pca, x='PC1', y='PC2', color='Race',
                               title="Analyse PCA - Structure génétique des populations")
                st.plotly_chart(fig, use_container_width=True)
                
                st.write(f"**Variance expliquée:**")
                for i, var in enumerate(pca_results['explained_variance'], 1):
                    st.write(f"PC{i}: {var*100:.1f}%")
    
    # Tab 4: SNP & QTL
    with genetic_tabs[3]:
        st.markdown("### 🧮 Analyse des SNP et QTL")
        
        snp_tabs = st.tabs(["SNP Database", "QTL Browser", "Association Studies", "Haplotypes"])
        
        with snp_tabs[0]:
            st.markdown("#### Base de données des Marqueurs SNP")
            
            # Recherche de SNP
            snp_search = st.text_input("Rechercher SNP (rsID ou position):", 
                                      placeholder="Ex: rs123456 ou chr5:123456")
            
            if snp_search:
                # Données simulées
                snp_data = {
                    'rs123456': {
                        'rs_id': 'rs123456',
                        'chromosome': '5',
                        'position': 123456,
                        'alleles': 'A/G',
                        'maf': 0.42,
                        'gene': 'MSTN',
                        'function': 'Missense',
                        'impact': 'Moderate',
                        'associated_trait': 'Muscle development'
                    },
                    'chr5:123456': {
                        'rs_id': 'rs123456',
                        'chromosome': '5',
                        'position': 123456,
                        'alleles': 'A/G',
                        'maf': 0.42,
                        'gene': 'MSTN',
                        'function': 'Missense',
                        'impact': 'Moderate',
                        'associated_trait': 'Muscle development'
                    }
                }
                
                result = snp_data.get(snp_search, {})
                
                if result:
                    st.markdown(f"##### SNP: **{result['rs_id']}**")
                    
                    col_snp1, col_snp2 = st.columns(2)
                    
                    with col_snp1:
                        st.write(f"**Chromosome:** {result['chromosome']}")
                        st.write(f"**Position:** {result['position']:,}")
                        st.write(f"**Allèles:** {result['alleles']}")
                        st.write(f"**MAF:** {result['maf']}")
                    
                    with col_snp2:
                        st.write(f"**Gène:** {result['gene']}")
                        st.write(f"**Fonction:** {result['function']}")
                        st.write(f"**Impact:** {result['impact']}")
                        st.write(f"**Trait associé:** {result['associated_trait']}")
                    
                    # Fréquence allélique par race
                    st.markdown("#### Fréquences par race")
                    
                    freq_data = pd.DataFrame({
                        'Race': ['Ouled_Djellal', 'Razè', 'Hamra', 'Dman', 'Saharienne'],
                        'Fréquence A': [0.8, 0.6, 0.7, 0.9, 0.5],
                        'Fréquence G': [0.2, 0.4, 0.3, 0.1, 0.5]
                    })
                    
                    fig = px.bar(freq_data, x='Race', y=['Fréquence A', 'Fréquence G'],
                                title="Fréquences alléliques par race",
                                barmode='group')
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("SNP non trouvé dans la base de données")
        
        with snp_tabs[1]:
            st.markdown("#### Navigateur QTL Ovin")
            
            # Liste des QTL connus
            qtl_db = [
                {'nom': 'QTL_MSTN', 'chromosome': '2', 'trait': 'Masse musculaire', 'lod': 12.5, 'var': 0.15},
                {'nom': 'QTL_MILK', 'chromosome': '6', 'trait': 'Production laitière', 'lod': 8.2, 'var': 0.08},
                {'nom': 'QTL_FAT', 'chromosome': '14', 'trait': 'Matière grasse', 'lod': 9.7, 'var': 0.12},
                {'nom': 'QTL_LITTER', 'chromosome': '10', 'trait': 'Taille portée', 'lod': 7.3, 'var': 0.07},
            ]
            
            df_qtl = pd.DataFrame(qtl_db)
            st.dataframe(df_qtl)
            
            # Graphique des LOD scores
            fig = px.bar(df_qtl, x='nom', y='lod', color='trait',
                        title="LOD Scores des QTL ovins")
            st.plotly_chart(fig, use_container_width=True)
    
    # Tab 5: Pedigrees
    with genetic_tabs[4]:
        st.markdown("### 🌳 Analyse de Pedigrees")
        
        pedigree_tabs = st.tabs(["Entrée", "Visualisation", "Calculs", "Optimisation"])
        
        with pedigree_tabs[0]:
            st.markdown("#### Saisie du Pedigree")
            
            pedigree_input = st.text_area(
                "Entrez le pedigree (Animal, Père, Mère):",
                height=200,
                placeholder="Animal1, Père1, Mère1\nAnimal2, Père1, Mère2\nAnimal3, Père2, Mère3\n..."
            )
            
            if st.button("📊 Analyser le pedigree"):
                if pedigree_input:
                    # Parser le pedigree
                    pedigree_data = []
                    for line in pedigree_input.strip().split('\n'):
                        parts = [p.strip() for p in line.split(',')]
                        if len(parts) >= 3:
                            pedigree_data.append(tuple(parts[:3]))
                    
                    if pedigree_data:
                        analysis = GeneticAnalyzer.pedigree_analyzer(pedigree_data)
                        
                        col_ped1, col_ped2, col_ped3 = st.columns(3)
                        
                        with col_ped1:
                            st.metric("Animaux", analysis['total_animals'])
                        with col_ped2:
                            st.metric("Relations", analysis['total_relations'])
                        with col_ped3:
                            st.metric("Consanguinité moyenne", f"{analysis['average_inbreeding']:.4f}")
                        
                        # Table des coefficients
                        df_inbreeding = pd.DataFrame(
                            list(analysis['inbreeding_coefficients'].items()),
                            columns=['Animal', 'Coefficient']
                        )
                        st.dataframe(df_inbreeding.sort_values('Coefficient', ascending=False))
        
        with pedigree_tabs[1]:
            st.markdown("#### Visualisation du Pedigree")
            st.info("Visualisation graphique en cours de développement...")
            
            # Graphique simplifié
            if 'pedigree_data' in locals() and pedigree_data:
                st.write("**Structure du pedigree:**")
                
                # Créer un graphe simple
                G = nx.DiGraph()
                for animal, sire, dam in pedigree_data:
                    if sire:
                        G.add_edge(sire, animal)
                    if dam:
                        G.add_edge(dam, animal)
                
                # Calculer les positions pour visualisation
                pos = nx.spring_layout(G, seed=42)
                
                fig, ax = plt.subplots(figsize=(10, 8))
                nx.draw(G, pos, with_labels=True, ax=ax, node_size=500, 
                       node_color='lightblue', font_size=8)
                st.pyplot(fig)
    
    # Tab 6: GWAS
    with genetic_tabs[5]:
        st.markdown("### 📈 Genome-Wide Association Study (GWAS)")
        
        st.markdown("""
        #### Analyse GWAS Complète
        
        Cette section permet d'effectuer des analyses GWAS pour identifier
        des marqueurs génétiques associés à des traits d'intérêt.
        
        **Traits disponibles:**
        - Production laitière
        - Teneur en matière grasse
        - Croissance musculaire
        - Fertilité
        - Résistance aux maladies
        """)
        
        trait_selected = st.selectbox("Trait à étudier:", 
                                     ["Production laitière", "Matière grasse", "Masse musculaire", "Fertilité"])
        
        if st.button("🎯 Lancer l'analyse GWAS"):
            with st.spinner("Analyse GWAS en cours... Cela peut prendre quelques minutes"):
                # Simulation d'analyse GWAS
                np.random.seed(42)
                
                # Générer des données simulées
                n_snps = 1000
                chromosomes = np.random.choice(range(1, 27), n_snps)
                positions = np.random.randint(1, 1000000, n_snps)
                p_values = np.random.exponential(0.1, n_snps)
                p_values = np.clip(p_values, 0, 1)
                
                # QQ-plot
                st.markdown("#### QQ-Plot")
                
                fig, ax = plt.subplots(figsize=(8, 6))
                stats.probplot(p_values, dist="uniform", plot=ax)
                ax.set_title("QQ-Plot - Distribution des p-values")
                st.pyplot(fig)
                
                # Manhattan plot
                st.markdown("#### Manhattan Plot")
                
                df_gwas = pd.DataFrame({
                    'CHR': chromosomes,
                    'POS': positions,
                    'P': p_values
                })
                
                # Créer le Manhattan plot
                df_gwas['-log10(P)'] = -np.log10(df_gwas['P'])
                
                fig = px.scatter(df_gwas, x='POS', y='-log10(P)', color='CHR',
                               title=f"Manhattan Plot - GWAS pour {trait_selected}",
                               labels={'POS': 'Position', '-log10(P)': '-log₁₀(p-value)'})
                
                # Ligne de significativité
                fig.add_hline(y=-np.log10(0.05/n_snps), line_dash="dash", 
                            line_color="red", annotation_text="Seuil de Bonferroni")
                fig.add_hline(y=-np.log10(0.001), line_dash="dot", 
                            line_color="orange", annotation_text="p < 0.001")
                
                st.plotly_chart(fig, use_container_width=True)
                
                # SNPs significatifs
                threshold = 0.05 / n_snps
                significant_snps = df_gwas[df_gwas['P'] < threshold]
                
                st.metric("SNPs significatifs (Bonferroni)", len(significant_snps))
    
    # Tab 7: Import/Export
    with genetic_tabs[6]:
        st.markdown("### 💾 Import/Export de Données Génétiques")
        
        format_tabs = st.tabs(["VCF", "PLINK", "FASTA", "GenBank"])
        
        with format_tabs[0]:
            st.markdown("#### Format VCF (Variant Call Format)")
            
            uploaded_vcf = st.file_uploader("Importer un fichier VCF", type=['vcf', 'vcf.gz'])
            
            if uploaded_vcf:
                content = uploaded_vcf.getvalue().decode()[:1000]
                st.success(f"Fichier VCF chargé: {uploaded_vcf.name}")
                st.code(content, language='vcf')
                
                if st.button("Analyser le VCF"):
                    st.info("Analyse VCF en cours de développement...")
            
            # Export VCF
            st.markdown("#### Exporter en VCF")
            if st.button("📥 Générer un VCF d'exemple"):
                # Créer un VCF exemple
                vcf_header = """##fileformat=VCFv4.2
##source=OvinManagerPro
##reference=Ovis_aries_3.1
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
"""
                
                vcf_data = """2\t123456\trs123456\tA\tG\t100\tPASS\tAF=0.42;DP=50\tGT\t0/0\t0/1\t1/1
6\t654321\trs654321\tC\tT\t150\tPASS\tAF=0.18;DP=60\tGT\t0/1\t0/0\t1/1
14\t987654\trs987654\tG\tA\t200\tPASS\tAF=0.33;DP=45\tGT\t1/1\t0/1\t0/0
"""
                
                vcf_content = vcf_header + vcf_data
                
                st.download_button(
                    label="📥 Télécharger VCF exemple",
                    data=vcf_content,
                    file_name="ovin_genotypes_example.vcf",
                    mime="text/plain"
                )
        
        with format_tabs[1]:
            st.markdown("#### Format PLINK")
            st.info("Support PLINK (.bed/.bim/.fam) en développement...")
        
        with format_tabs[2]:
            st.markdown("#### Format FASTA")
            
            fasta_input = st.text_area("Entrez des séquences FASTA:", height=200)
            
            if fasta_input and st.button("Convertir FASTA"):
                st.code(fasta_input, language='fasta')
        
        with format_tabs[3]:
            st.markdown("#### Format GenBank")
            st.info("Import/Export GenBank en développement...")

# ========== NAVIGATION PRINCIPALE ==========

# Titre principal
st.markdown('<h1 class="main-header">🐑 Ovin Manager Pro - Génétique Avancée</h1>', unsafe_allow_html=True)
st.markdown("""
*Application scientifique complète de gestion ovine avec module génétique professionnel*
""")

# Sidebar - Navigation
with st.sidebar:
    st.markdown("### 📍 Navigation Génétique")
    
    page = st.radio(
        "Menu Principal",
        ["🏠 Tableau de Bord", 
         "📊 Gestion des Brebis", 
         "🧬 Génétique Avancée",      # Module génétique développé
         "🎯 Scoring Phénotypique", 
         "📝 Formulaires Standardisés",
         "🥛 Analyse Lait",
         "📐 Morphométrie 3D",
         "🤰 Suivi Gestation", 
         "📈 Statistiques Avancées",
         "⚙️ Paramètres"]
    )
    
    st.markdown("---")
    st.markdown("### 📊 Métriques Génétiques")
    
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM donnees_genomiques")
    seq_count = cursor.fetchone()[0]
    st.metric("Séquences stockées", seq_count)
    
    cursor.execute("SELECT COUNT(DISTINCT gene_nom) FROM donnees_genomiques")
    gene_count = cursor.fetchone()[0]
    st.metric("Gènes uniques", gene_count)
    
    cursor.execute("SELECT COUNT(*) FROM snp_marqueurs")
    snp_count = cursor.fetchone()[0]
    st.metric("Marqueurs SNP", snp_count)

# Navigation principale
if page == "🏠 Tableau de Bord":
    # (Code du tableau de bord existant à ajouter ici)
    st.info("Tableau de bord - À intégrer")
elif page == "📊 Gestion des Brebis":
    # (Code gestion brebis existant)
    st.info("Gestion des brebis - À intégrer")
elif page == "🧬 Génétique Avancée":
    afficher_genetique_avancee()
elif page == "🎯 Scoring Phénotypique":
    # (Code scoring phénotypique existant)
    st.info("Scoring phénotypique - À intégrer")
elif page == "📝 Formulaires Standardisés":
    # (Code formulaires existant)
    st.info("Formulaires standardisés - À intégrer")
elif page == "🥛 Analyse Lait":
    # (Code analyse lait existant)
    st.info("Analyse laitière - À intégrer")
elif page == "📐 Morphométrie 3D":
    # (Code morphométrie existant)
    st.info("Morphométrie 3D - À intégrer")
elif page == "🤰 Suivi Gestation":
    # (Code gestation existant)
    st.info("Suivi gestation - À intégrer")
elif page == "📈 Statistiques Avancées":
    # (Code statistiques existant)
    st.info("Statistiques avancées - À intégrer")
elif page == "⚙️ Paramètres":
    # (Code paramètres existant)
    st.info("Paramètres - À intégrer")

# Pied de page
st.markdown("---")
st.caption("🧬 Ovin Manager Pro v3.0 - Module Génétique Avancé | © 2024")
