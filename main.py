"""
OVIN MANAGER PRO - Version Professionnelle avec Génétique Avancée
Sans dépendances problématiques pour Streamlit Cloud
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
import math
import random
from collections import Counter

# ========== CONFIGURATION ==========

st.set_page_config(
    page_title="Ovin Manager Pro - Génétique Avancée",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS personnalisé avec animations
st.markdown("""
<style>
    .main-header {
        font-size: 2.8rem;
        color: #2E7D32;
        text-align: center;
        margin-bottom: 1rem;
        background: linear-gradient(90deg, #2E7D32, #4CAF50);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .section-header {
        font-size: 2rem;
        color: #388E3C;
        margin-top: 2rem;
        margin-bottom: 1rem;
        padding-bottom: 10px;
        border-bottom: 3px solid #4CAF50;
    }
    .gene-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        border-radius: 15px;
        margin: 15px 0;
        box-shadow: 0 10px 20px rgba(0,0,0,0.1);
    }
    .dna-sequence {
        font-family: 'Courier New', monospace;
        background: #1a1a2e;
        color: #00ff88;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
        font-size: 1em;
        letter-spacing: 2px;
        overflow-x: auto;
        border-left: 5px solid #00ff88;
    }
    .snp-card {
        background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        color: white;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
    }
    .metric-enhanced {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        border-left: 5px solid #4CAF50;
    }
</style>
""", unsafe_allow_html=True)

# ========== BASE DE DONNÉES AVANCÉE ==========

def init_advanced_db():
    """Initialise une base de données génétique avancée"""
    conn = sqlite3.connect('ovin_genetic.db', check_same_thread=False)
    cursor = conn.cursor()
    
    # Table brebis détaillée
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ovins (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant TEXT UNIQUE NOT NULL,
            nom TEXT,
            race TEXT,
            sous_race TEXT,
            sexe TEXT CHECK(sexe IN ('F', 'M')),
            date_naissance DATE,
            poids FLOAT,
            score_condition INTEGER CHECK(score_condition BETWEEN 1 AND 5),
            couleur_robe TEXT,
            cornes BOOLEAN,
            notes TEXT,
            mere_id TEXT,
            pere_id TEXT,
            coefficient_consanguinite FLOAT DEFAULT 0.0,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Table génome
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS genome_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            ovins_id INTEGER,
            chromosome TEXT,
            position_start INTEGER,
            position_end INTEGER,
            gene_nom TEXT,
            gene_symbol TEXT,
            sequence_adn TEXT,
            genotype TEXT,
            allele1 TEXT,
            allele2 TEXT,
            variant_type TEXT,
            impact TEXT CHECK(impact IN ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')),
            quality_score FLOAT,
            read_depth INTEGER,
            allele_frequency FLOAT,
            date_analyse DATE,
            FOREIGN KEY (ovins_id) REFERENCES ovins(id)
        )
    ''')
    
    # Table marqueurs génétiques
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS genetic_markers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            rs_id TEXT UNIQUE,
            chromosome TEXT,
            position INTEGER,
            ref_allele TEXT,
            alt_allele TEXT,
            gene_symbol TEXT,
            gene_name TEXT,
            consequence TEXT,
            maaf_global FLOAT,
            maaf_ouled_djellal FLOAT,
            maaf_raze FLOAT,
            maaf_hamra FLOAT,
            heritabilite FLOAT,
            trait_associe TEXT,
            publication TEXT,
            date_ajout DATE
        )
    ''')
    
    # Table QTL (Quantitative Trait Loci)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS qtl_ovins (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nom TEXT,
            chromosome TEXT,
            start_position INTEGER,
            end_position INTEGER,
            trait TEXT,
            lod_score FLOAT,
            variance_explained FLOAT,
            peak_marker TEXT,
            genes_candidats TEXT,
            race_studiee TEXT,
            publication TEXT
        )
    ''')
    
    # Table analyses GWAS
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS gwas_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            trait TEXT,
            chromosome TEXT,
            position INTEGER,
            p_value FLOAT,
            effect_size FLOAT,
            marker_id TEXT,
            allele_effectif TEXT,
            beta FLOAT,
            se FLOAT,
            bonferroni_threshold FLOAT,
            significant BOOLEAN,
            date_analyse DATE
        )
    ''')
    
    # Table pedigrees
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pedigrees (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            animal_id TEXT UNIQUE,
            sire_id TEXT,
            dam_id TEXT,
            sire_sire TEXT,
            sire_dam TEXT,
            dam_sire TEXT,
            dam_dam TEXT,
            generation INTEGER,
            inbreeding_coefficient FLOAT,
            additive_relationship FLOAT
        )
    ''')
    
    # Peupler avec des données de démonstration
    peupler_donnees_demo(cursor)
    
    conn.commit()
    return conn

def peupler_donnees_demo(cursor):
    """Peuple la base avec des données de démonstration"""
    
    # Insérer des marqueurs génétiques de démonstration
    markers = [
        ('rs123456789', '2', 123456, 'A', 'G', 'MSTN', 'Myostatine', 'missense_variant', 0.42, 0.38, 0.45, 0.41, 0.85, 'Muscle mass', 'PMID:12345678'),
        ('rs987654321', '6', 654321, 'C', 'T', 'PRNP', 'Prion protein', 'synonymous_variant', 0.18, 0.15, 0.20, 0.19, 0.92, 'Scrapie resistance', 'PMID:87654321'),
        ('rs456789123', '14', 789123, 'G', 'A', 'DGAT1', 'Diacylglycerol acyltransferase', 'missense_variant', 0.33, 0.28, 0.35, 0.30, 0.45, 'Milk fat content', 'PMID:34567891'),
        ('rs321654987', '5', 321654, 'T', 'C', 'GDF9', 'Growth differentiation factor 9', 'frameshift_variant', 0.25, 0.22, 0.27, 0.24, 0.65, 'Litter size', 'PMID:21654987'),
        ('rs741852963', '10', 741852, 'A', 'T', 'BMPR1B', 'Bone morphogenetic protein receptor', 'stop_gained', 0.12, 0.10, 0.15, 0.11, 0.75, 'Fecundity', 'PMID:41852963')
    ]
    
    for marker in markers:
        try:
            cursor.execute('''
                INSERT OR IGNORE INTO genetic_markers 
                (rs_id, chromosome, position, ref_allele, alt_allele, gene_symbol, gene_name, 
                 consequence, maaf_global, maaf_ouled_djellal, maaf_raze, maaf_hamra, 
                 heritabilite, trait_associe, publication, date_ajout)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (*marker, date.today().isoformat()))
        except:
            pass
    
    # Insérer des QTL de démonstration
    qtls = [
        ('QTL_MSTN_1', '2', 1200000, 1300000, 'Muscle mass', 12.5, 0.15, 'rs123456789', 'MSTN, MYF5, MYOD1', 'Ouled Djellal', 'PMID:12345678'),
        ('QTL_MILK_1', '6', 650000, 700000, 'Milk yield', 8.2, 0.08, 'rs987654321', 'PRNP, CSN2, CSN3', 'Razè', 'PMID:87654321'),
        ('QTL_FAT_1', '14', 780000, 800000, 'Milk fat', 9.7, 0.12, 'rs456789123', 'DGAT1, FASN, SCD', 'Hamra', 'PMID:34567891'),
        ('QTL_FEC_1', '5', 320000, 330000, 'Fecundity', 7.3, 0.07, 'rs321654987', 'GDF9, BMP15, BMPR1B', 'Dman', 'PMID:21654987')
    ]
    
    for qtl in qtls:
        try:
            cursor.execute('''
                INSERT OR IGNORE INTO qtl_ovins 
                (nom, chromosome, start_position, end_position, trait, lod_score, 
                 variance_explained, peak_marker, genes_candidats, race_studiee, publication)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', qtl)
        except:
            pass

# Initialisation de la base de données
conn = init_advanced_db()

# ========== CLASSES GÉNÉTIQUES PERSONNALISÉES ==========

class GeneticAnalyzerCustom:
    """Analyseur génétique personnalisé sans dépendances externes"""
    
    @staticmethod
    def analyze_dna_sequence(sequence):
        """Analyse approfondie d'une séquence ADN"""
        seq = sequence.upper().replace(' ', '').replace('\n', '')
        
        if len(seq) == 0:
            return {"error": "Séquence vide"}
        
        # Composition détaillée
        bases = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, 'Other': 0}
        for base in seq:
            if base in bases:
                bases[base] += 1
            else:
                bases['Other'] += 1
        
        total = len(seq)
        
        # Calculs avancés
        gc_content = (bases['G'] + bases['C']) / total * 100
        at_content = 100 - gc_content
        gc_skew = (bases['G'] - bases['C']) / (bases['G'] + bases['C']) if (bases['G'] + bases['C']) > 0 else 0
        at_skew = (bases['A'] - bases['T']) / (bases['A'] + bases['T']) if (bases['A'] + bases['T']) > 0 else 0
        
        # Recherche de motifs
        motifs = {
            'start_codons': seq.count('ATG'),
            'stop_codons': seq.count('TAA') + seq.count('TAG') + seq.count('TGA'),
            'cpg_islands': GeneticAnalyzerCustom._find_cpg_islands(seq),
            'repeats': GeneticAnalyzerCustom._find_repeats(seq),
            'restriction_sites': GeneticAnalyzerCustom._find_restriction_sites(seq)
        }
        
        # Prédictions
        predictions = {
            'coding_potential': GeneticAnalyzerCustom._predict_coding_potential(seq),
            'secondary_structure': GeneticAnalyzerCustom._predict_secondary_structure(seq),
            'melting_temperature': GeneticAnalyzerCustom._calculate_melting_temp(seq),
            'molecular_weight': GeneticAnalyzerCustom._calculate_molecular_weight(seq)
        }
        
        return {
            'length': total,
            'composition': bases,
            'percentages': {
                'GC': round(gc_content, 2),
                'AT': round(at_content, 2),
                'GC_skew': round(gc_skew, 3),
                'AT_skew': round(at_skew, 3)
            },
            'motifs': motifs,
            'predictions': predictions,
            'checksum': hash(seq) % 1000000
        }
    
    @staticmethod
    def _find_cpg_islands(seq, window=200, threshold=0.5):
        """Détecte les îlots CpG"""
        islands = []
        for i in range(0, len(seq) - window + 1, window//2):
            subseq = seq[i:i+window]
            cpg_count = subseq.count('CG')
            gc_content = (subseq.count('G') + subseq.count('C')) / len(subseq)
            
            if cpg_count > 0 and gc_content > threshold:
                islands.append({
                    'start': i,
                    'end': i + window,
                    'cpg_density': cpg_count / window,
                    'gc_content': gc_content
                })
        return islands[:10]
    
    @staticmethod
    def _find_repeats(seq):
        """Trouve les séquences répétées"""
        repeats = {}
        
        # Recherche de microsatellites
        for repeat_len in [2, 3, 4]:
            max_repeats = 0
            best_motif = ''
            
            for i in range(len(seq) - repeat_len * 3 + 1):
                motif = seq[i:i+repeat_len]
                count = 1
                
                # Compter les répétitions consécutives
                j = i + repeat_len
                while j <= len(seq) - repeat_len and seq[j:j+repeat_len] == motif:
                    count += 1
                    j += repeat_len
                
                if count > max_repeats and count >= 3:
                    max_repeats = count
                    best_motif = motif
            
            if max_repeats > 0:
                repeats[f'{repeat_len}-mer'] = {
                    'motif': best_motif,
                    'repeats': max_repeats,
                    'length': len(best_motif) * max_repeats
                }
        
        return repeats
    
    @staticmethod
    def _find_restriction_sites(seq):
        """Trouve les sites de restriction courants"""
        enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XbaI': 'TCTAGA',
            'NotI': 'GCGGCCGC',
            'SacI': 'GAGCTC',
            'PstI': 'CTGCAG'
        }
        
        sites = []
        for enzyme, site in enzymes.items():
            positions = []
            for i in range(len(seq) - len(site) + 1):
                if seq[i:i+len(site)] == site:
                    positions.append(i)
            
            if positions:
                sites.append({
                    'enzyme': enzyme,
                    'recognition_site': site,
                    'positions': positions,
                    'count': len(positions)
                })
        
        return sites
    
    @staticmethod
    def _predict_coding_potential(seq):
        """Prédit le potentiel de codage"""
        if len(seq) < 100:
            return {'score': 0, 'confidence': 'low', 'likely_coding': False}
        
        # Score basé sur plusieurs facteurs
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        
        # Vérifier les cadres de lecture
        orf_scores = []
        for frame in range(3):
            orf_length = 0
            max_orf = 0
            
            for i in range(frame, len(seq)-2, 3):
                codon = seq[i:i+3]
                if codon == 'ATG':  # Start codon
                    orf_length = 3
                elif codon in ['TAA', 'TAG', 'TGA']:  # Stop codon
                    if orf_length > max_orf:
                        max_orf = orf_length
                    orf_length = 0
                else:
                    orf_length += 3
            
            orf_scores.append(max_orf)
        
        max_orf_score = max(orf_scores) / len(seq)
        
        # Score composite
        coding_score = (gc_content * 0.3 + max_orf_score * 0.7) * 100
        
        return {
            'score': round(coding_score, 2),
            'confidence': 'high' if coding_score > 70 else 'medium' if coding_score > 40 else 'low',
            'likely_coding': coding_score > 50,
            'max_orf_length': max(orf_scores),
            'best_frame': orf_scores.index(max(orf_scores)) + 1
        }
    
    @staticmethod
    def _predict_secondary_structure(seq):
        """Prédit la structure secondaire simplifiée"""
        gc_content = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
        
        if gc_content > 0.6:
            structure = "GC-rich, likely stable with strong secondary structure"
            stability = "high"
        elif gc_content > 0.4:
            structure = "Balanced GC/AT, moderate structure"
            stability = "medium"
        else:
            structure = "AT-rich, less structured, more flexible"
            stability = "low"
        
        # Calcul de la température de fusion approximative
        tm = 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))
        
        return {
            'prediction': structure,
            'stability': stability,
            'gc_content': round(gc_content * 100, 1),
            'melting_temperature': tm,
            'stem_loop_potential': round(gc_content * 100, 1)
        }
    
    @staticmethod
    def _calculate_melting_temp(seq):
        """Calcule la température de fusion (formule simplifiée)"""
        return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))
    
    @staticmethod
    def _calculate_molecular_weight(seq):
        """Calcule le poids moléculaire"""
        weights = {'A': 313.21, 'T': 304.20, 'C': 289.18, 'G': 329.21, 'N': 300.00}
        total = sum(weights.get(base, 300.00) for base in seq)
        return round(total / 1000, 2)  # en kDa
    
    @staticmethod
    def translate_sequence(seq):
        """Traduit une séquence ADN en protéine"""
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
        }
        
        protein = ""
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i+3]
            protein += codon_table.get(codon, 'X')
        
        return protein
    
    @staticmethod
    def calculate_hardy_weinberg(genotypes):
        """Calcule l'équilibre de Hardy-Weinberg"""
        if not genotypes:
            return None
        
        counts = Counter(genotypes)
        total = sum(counts.values())
        
        # Extraire les allèles
        alleles = []
        for genotype in genotypes:
            alleles.extend([genotype[0], genotype[1]])
        
        allele_counts = Counter(alleles)
        allele_freq = {allele: count/(total*2) for allele, count in allele_counts.items()}
        
        # Calcul des fréquences attendues
        expected_counts = {}
        alleles_list = list(allele_freq.keys())
        
        for i, a1 in enumerate(alleles_list):
            for a2 in alleles_list[i:]:
                genotype = ''.join(sorted([a1, a2]))
                freq = allele_freq[a1] * allele_freq[a2]
                if a1 != a2:
                    freq *= 2
                expected_counts[genotype] = round(freq * total, 2)
        
        # Test du chi²
        chi2 = 0
        all_genotypes = set(list(counts.keys()) + list(expected_counts.keys()))
        
        for genotype in all_genotypes:
            observed = counts.get(genotype, 0)
            expected = expected_counts.get(genotype, 0)
            if expected > 0:
                chi2 += ((observed - expected) ** 2) / expected
        
        # Degrés de liberté
        df = len(alleles_list) * (len(alleles_list) - 1) // 2
        
        return {
            'observed': dict(counts),
            'expected': expected_counts,
            'allele_frequencies': allele_freq,
            'chi_squared': round(chi2, 4),
            'degrees_of_freedom': df,
            'p_value': round(math.exp(-chi2/2) * (chi2**(df/2-1)) / (2**(df/2) * math.gamma(df/2)), 6),
            'in_equilibrium': chi2 < 3.841  # p > 0.05 pour df=1
        }

class PopulationGeneticsCustom:
    """Génétique des populations personnalisée"""
    
    @staticmethod
    def calculate_f_statistics(subpopulations):
        """Calcule les F-statistiques (Fis, Fit, Fst)"""
        # Sous-populations: liste de listes de génotypes
        if not subpopulations or len(subpopulations) < 2:
            return None
        
        # Calcul global
        all_genotypes = []
        for pop in subpopulations:
            all_genotypes.extend(pop)
        
        # Fréquences alléliques globales
        all_alleles = []
        for genotype in all_genotypes:
            all_alleles.extend([genotype[0], genotype[1]])
        
        global_allele_freq = dict(Counter(all_alleles))
        total_alleles = len(all_alleles)
        for allele in global_allele_freq:
            global_allele_freq[allele] /= total_alleles
        
        # Calcul par population
        pop_stats = []
        for pop in subpopulations:
            pop_alleles = []
            for genotype in pop:
                pop_alleles.extend([genotype[0], genotype[1]])
            
            pop_allele_freq = dict(Counter(pop_alleles))
            pop_total = len(pop_alleles)
            for allele in pop_allele_freq:
                pop_allele_freq[allele] = pop_allele_freq.get(allele, 0) / pop_total
            
            # Hétérozygotie observée
            het_obs = sum(1 for g in pop if g[0] != g[1]) / len(pop) if len(pop) > 0 else 0
            
            # Hétérozygotie attendue
            het_exp = 1 - sum(f**2 for f in pop_allele_freq.values())
            
            pop_stats.append({
                'sample_size': len(pop),
                'allele_frequencies': pop_allele_freq,
                'heterozygosity_observed': het_obs,
                'heterozygosity_expected': het_exp
            })
        
        # Calcul F-statistiques
        Hs = np.mean([p['heterozygosity_observed'] for p in pop_stats])
        Ht = 1 - sum(f**2 for f in global_allele_freq.values())
        
        Fis = 1 - (Hs / Ht) if Ht > 0 else 0
        Fst = (Ht - Hs) / Ht if Ht > 0 else 0
        Fit = 1 - (Hs / Ht) if Ht > 0 else 0
        
        return {
            'subpopulations': len(subpopulations),
            'total_individuals': len(all_genotypes),
            'fis_inbreeding': round(Fis, 4),
            'fst_population_differentiation': round(Fst, 4),
            'fit_total_inbreeding': round(Fit, 4),
            'heterozygosity_within': round(Hs, 4),
            'heterozygosity_total': round(Ht, 4),
            'population_stats': pop_stats
        }

# ========== PAGES PROFESSIONNELLES ==========

def page_genetique_pro():
    """Page génétique professionnelle"""
    st.markdown('<h1 class="main-header">🧬 LABORATOIRE DE GÉNÉTIQUE OVINE</h1>', unsafe_allow_html=True)
    st.markdown("*Outils professionnels d'analyse génomique pour la recherche et la sélection*")
    
    # Menu génétique avancé
    genetique_tabs = st.tabs([
        "🧬 ANALYSE GÉNOME", 
        "🔍 RECHERCHE GÈNES", 
        "📊 GWAS & QTL", 
        "🧮 POP. GENETICS",
        "🌳 PEDIGREES",
        "💾 IMPORT/EXPORT"
    ])
    
    # Tab 1: Analyse Génome
    with genetique_tabs[0]:
        st.markdown('<h2 class="section-header">🧬 Analyse Complète de Séquences Génomiques</h2>', unsafe_allow_html=True)
        
        col_seq1, col_seq2 = st.columns([3, 1])
        
        with col_seq1:
            sequence_input = st.text_area(
                "**ENTREZ VOTRE SÉQUENCE ADN** (FASTA ou format brut):",
                height=300,
                placeholder=">Sequence_Ovis_aries_gene_MSTN\nATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGCCTGGGCCGCATGGCCAAAGAGATCAAG...\n\nOu directement la séquence:\nATCGATCGATCG..."
            )
            
            # Options d'analyse
            analysis_options = st.multiselect(
                "**SÉLECTIONNEZ LES ANALYSES:**",
                [
                    "Composition nucléotidique détaillée",
                    "Recherche de motifs fonctionnels",
                    "Prédiction de structure secondaire",
                    "Sites de restriction",
                    "Potentiel de codage (ORF)",
                    "Traduction en protéine",
                    "Analyse de conservation"
                ],
                default=["Composition nucléotidique détaillée", "Recherche de motifs fonctionnels"]
            )
        
        with col_seq2:
            st.markdown("""
            <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        color: white; padding: 20px; border-radius: 15px; margin-bottom: 20px;'>
                <h4>📚 RESSOURCES</h4>
                <p>• Référence: Ovis aries v4.0</p>
                <p>• 27 chromosomes</p>
                <p>• 2.7 Gb de génome</p>
                <p>• ~20,000 gènes annotés</p>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            # Exemples de séquences
            st.markdown("**EXEMPLES:**")
            examples = {
                "MSTN (Myostatine)": "ATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGCCTGGGCCGCATGGCCAAAGAGATCAAG",
                "PRNP (Prion)": "ATGCGAACCTTGGAGGCGGTGGCTTCCTCGCTGCTGGTAGCGGCGGTGGCGGTGGCTTCCTCGCTGGTGGTAGC",
                "DGAT1": "ATGGAGAGCGCCGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAG"
            }
            
            for gene, seq in examples.items():
                if st.button(f"📥 {gene}"):
                    st.session_state.current_sequence = seq
                    st.rerun()
        
        # Bouton d'analyse
        if st.button("🚀 LANCER L'ANALYSE GÉNOMIQUE", type="primary", use_container_width=True):
            if sequence_input:
                with st.spinner("🔬 Analyse en cours... Cela peut prendre quelques secondes"):
                    # Nettoyer la séquence
                    lines = sequence_input.strip().split('\n')
                    sequence = ''.join([line for line in lines if not line.startswith('>')])
                    sequence = sequence.upper().replace(' ', '').replace('\n', '')
                    
                    if len(sequence) > 0:
                        # Analyse complète
                        results = GeneticAnalyzerCustom.analyze_dna_sequence(sequence)
                        
                        if "error" not in results:
                            st.success(f"✅ ANALYSE TERMINÉE ! Séquence de {results['length']:,} pb analysée")
                            
                            # Métriques principales
                            col_metrics1, col_metrics2, col_metrics3, col_metrics4 = st.columns(4)
                            
                            with col_metrics1:
                                st.markdown("""
                                <div class='metric-enhanced'>
                                    <h3>📏 LONGUEUR</h3>
                                    <h2>{:,} pb</h2>
                                </div>
                                """.format(results['length']), unsafe_allow_html=True)
                            
                            with col_metrics2:
                                st.markdown(f"""
                                <div class='metric-enhanced'>
                                    <h3>🧬 % GC</h3>
                                    <h2>{results['percentages']['GC']}%</h2>
                                </div>
                                """, unsafe_allow_html=True)
                            
                            with col_metrics3:
                                st.markdown(f"""
                                <div class='metric-enhanced'>
                                    <h3>🔬 TM</h3>
                                    <h2>{results['predictions']['melting_temperature']}°C</h2>
                                </div>
                                """, unsafe_allow_html=True)
                            
                            with col_metrics4:
                                st.markdown(f"""
                                <div class='metric-enhanced'>
                                    <h3>⚖️ POIDS</h3>
                                    <h2>{results['predictions']['molecular_weight']} kDa</h2>
                                </div>
                                """, unsafe_allow_html=True)
                            
                            # Graphiques et résultats détaillés
                            result_detail_tabs = st.tabs(["📊 COMPOSITION", "🎯 MOTIFS", "🧪 PRÉDICTIONS", "📋 RAPPORT"])
                            
                            with result_detail_tabs[0]:
                                # Graphique de composition
                                df_composition = pd.DataFrame({
                                    'Base': list(results['composition'].keys()),
                                    'Count': list(results['composition'].values()),
                                    'Percentage': [v/results['length']*100 for v in results['composition'].values()]
                                })
                                
                                fig1 = px.bar(df_composition, x='Base', y='Count', color='Base',
                                            title="Composition Nucléotidique Détaillée")
                                st.plotly_chart(fig1, use_container_width=True)
                                
                                # Diagramme circulaire
                                fig2 = px.pie(df_composition[df_composition['Base'] != 'Other'], 
                                            values='Count', names='Base', 
                                            title="Distribution des Bases")
                                st.plotly_chart(fig2, use_container_width=True)
                            
                            with result_detail_tabs[1]:
                                if results['motifs']:
                                    # Sites de restriction
                                    if results['motifs']['restriction_sites']:
                                        st.markdown("### 🧪 SITES DE RESTRICTION DÉTECTÉS")
                                        for site in results['motifs']['restriction_sites']:
                                            with st.expander(f"🔪 {site['enzyme']} ({site['recognition_site']})"):
                                                st.write(f"**Nombre:** {site['count']}")
                                                st.write(f"**Positions:** {site['positions'][:5]}{'...' if len(site['positions']) > 5 else ''}")
                                    
                                    # Répétitions
                                    if results['motifs']['repeats']:
                                        st.markdown("### 🔁 MICROSATELLITES")
                                        for repeat_type, repeat_info in results['motifs']['repeats'].items():
                                            st.write(f"**{repeat_type}:** {repeat_info['motif']} × {repeat_info['repeats']}")
                                    
                                    # Îlots CpG
                                    if results['motifs']['cpg_islands']:
                                        st.markdown("### 🏝️ ÎLOTS CpG")
                                        df_cpg = pd.DataFrame(results['motifs']['cpg_islands'])
                                        st.dataframe(df_cpg)
                            
                            with result_detail_tabs[2]:
                                # Prédiction de structure
                                st.markdown(f"""
                                <div class='gene-card'>
                                    <h3>🧬 PRÉDICTION DE STRUCTURE</h3>
                                    <p><strong>Type:</strong> {results['predictions']['secondary_structure']['prediction']}</p>
                                    <p><strong>Stabilité:</strong> {results['predictions']['secondary_structure']['stability'].upper()}</p>
                                    <p><strong>Température de fusion:</strong> {results['predictions']['secondary_structure']['melting_temperature']}°C</p>
                                </div>
                                """, unsafe_allow_html=True)
                                
                                # Potentiel de codage
                                coding = results['predictions']['coding_potential']
                                st.markdown(f"""
                                <div class='snp-card'>
                                    <h3>🧪 POTENTIEL DE CODAGE</h3>
                                    <p><strong>Score:</strong> {coding['score']}/100</p>
                                    <p><strong>Confiance:</strong> {coding['confidence'].upper()}</p>
                                    <p><strong>ORF le plus long:</strong> {coding['max_orf_length']} pb (cadre {coding['best_frame']})</p>
                                    <p><strong>Probablement codant:</strong> {'✅ OUI' if coding['likely_coding'] else '❌ NON'}</p>
                                </div>
                                """, unsafe_allow_html=True)
                                
                                # Traduction
                                if len(sequence) >= 3:
                                    protein = GeneticAnalyzerCustom.translate_sequence(sequence[:300])
                                    st.markdown("### 🧫 TRADUCTION PROTÉIQUE (100 premiers acides aminés)")
                                    st.markdown(f'<div class="dna-sequence">{protein[:100]}</div>', unsafe_allow_html=True)
                            
                            with result_detail_tabs[3]:
                                # Rapport complet
                                st.download_button(
                                    label="📥 TÉLÉCHARGER LE RAPPORT COMPLET (JSON)",
                                    data=json.dumps(results, indent=2),
                                    file_name=f"rapport_genomique_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                                    mime="application/json",
                                    type="primary"
                                )
                                
                                st.code(json.dumps(results, indent=2), language='json')
                        else:
                            st.error("❌ ERREUR: " + results["error"])
                    else:
                        st.warning("⚠️ Veuillez entrer une séquence valide")
            else:
                st.warning("⚠️ Veuillez entrer une séquence ADN")
    
    # Tab 2: Recherche de gènes
    with genetique_tabs[1]:
        st.markdown('<h2 class="section-header">🔍 BASE DE DONNÉES DES GÈNES OVINS</h2>', unsafe_allow_html=True)
        
        # Recherche de gènes
        col_search1, col_search2 = st.columns([3, 1])
        
        with col_search1:
            gene_query = st.text_input("**RECHERCHEZ UN GÈNE:**", placeholder="Ex: MSTN, PRNP, DGAT1, GDF9, BMPR1B...")
        
        with col_search2:
            search_type = st.selectbox("**TYPE DE RECHERCHE:**", ["Par nom", "Par fonction", "Par chromosome", "Par phénotype"])
        
        # Base de données des gènes ovins
        genes_database = {
            "MSTN": {
                "nom_complet": "Myostatine",
                "synonymes": ["GDF8", "Growth Differentiation Factor 8"],
                "chromosome": "2",
                "position": "chr2:6,254,871-6,265,123",
                "fonction": "Régulateur négatif de la croissance musculaire squelettique",
                "phénotype": ["Hypertrophie musculaire", "Double-muscling", "Accroissement de la masse maigre"],
                "mutations": [
                    {"id": "g.6723G>A", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "c.939G>A", "type": "SNP", "effet": "Synonyme", "impact": "Faible"},
                    {"id": "11-bp dél", "type": "Délétion", "effet": "Frameshift", "impact": "Élevé"}
                ],
                "heritabilite": 0.85,
                "qtl_associes": ["QTL_MSTN_1", "QTL_MUSCLE_2"],
                "publications": ["PMID:12345678", "PMID:23456789"],
                "sequence_exemple": "ATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGCCTGGGCCGCATGG"
            },
            "PRNP": {
                "nom_complet": "Protéine Prion",
                "synonymes": ["PrP", "CD230"],
                "chromosome": "13",
                "position": "chr13:42,316,543-42,328,976",
                "fonction": "Protéine membranaire, rôle dans la susceptibilité aux encéphalopathies spongiformes",
                "phénotype": ["Résistance à la tremblante", "Tolérance aux ESST", "Longévité"],
                "mutations": [
                    {"id": "codon 136", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "codon 154", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "codon 171", "type": "SNP", "effet": "Missense", "impact": "Élevé"}
                ],
                "heritabilite": 0.92,
                "qtl_associes": ["QTL_SCRAPIE_1"],
                "publications": ["PMID:34567891"],
                "sequence_exemple": "ATGCGAACCTTGGAGGCGGTGGCTTCCTCGCTGCTGGTAGCGGCGGTGGCGGTGGCTTCCTCGCT"
            },
            "DGAT1": {
                "nom_complet": "Diacylglycérol acyltransférase 1",
                "synonymes": ["ARGP1", "DGAT"],
                "chromosome": "14",
                "position": "chr14:21,894,765-21,912,345",
                "fonction": "Enzyme clé dans la biosynthèse des triglycérides",
                "phénotype": ["Teneur en matière grasse du lait", "Rendement fromager", "Énergie du lait"],
                "mutations": [
                    {"id": "K232A", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "c.10433C>T", "type": "SNP", "effet": "Synonyme", "impact": "Faible"}
                ],
                "heritabilite": 0.45,
                "qtl_associes": ["QTL_MILKFAT_1", "QTL_MILK_3"],
                "publications": ["PMID:45678912"],
                "sequence_exemple": "ATGGAGAGCGCCGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAG"
            },
            "GDF9": {
                "nom_complet": "Growth Differentiation Factor 9",
                "synonymes": ["GDF-9"],
                "chromosome": "5",
                "position": "chr5:88,456,123-88,462,456",
                "fonction": "Facteur de croissance impliqué dans la folliculogenèse et la fertilité femelle",
                "phénotype": ["Prolificité", "Taille de portée", "Fertilité"],
                "mutations": [
                    {"id": "G1", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "G4", "type": "SNP", "effet": "Missense", "impact": "Élevé"},
                    {"id": "G8", "type": "Indel", "effet": "Frameshift", "impact": "Élevé"}
                ],
                "heritabilite": 0.65,
                "qtl_associes": ["QTL_FECUNDITY_1", "QTL_LITTER_2"],
                "publications": ["PMID:56789123"],
                "sequence_exemple": "ATGGGCCCCGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAG"
            }
        }
        
        # Affichage des résultats
        if gene_query:
            gene_key = gene_query.upper()
            if gene_key in genes_database:
                gene_info = genes_database[gene_key]
                
                st.markdown(f"""
                <div class='gene-card'>
                    <h2>🧬 {gene_key} - {gene_info['nom_complet']}</h2>
                    <p><strong>📍 Localisation:</strong> {gene_info['chromosome']} ({gene_info['position']})</p>
                    <p><strong>🎯 Héritabilité:</strong> {gene_info['heritabilite']}</p>
                </div>
                """, unsafe_allow_html=True)
                
                # Onglets d'information
                gene_tabs = st.tabs(["📋 INFORMATIONS", "🧬 MUTATIONS", "📊 PHÉNOTYPES", "🔗 RESSOURCES"])
                
                with gene_tabs[0]:
                    col_gene1, col_gene2 = st.columns(2)
                    
                    with col_gene1:
                        st.markdown("### 📖 DESCRIPTION")
                        st.write(gene_info['fonction'])
                        
                        st.markdown("### 🔤 SYNONYMES")
                        st.write(", ".join(gene_info['synonymes']))
                    
                    with col_gene2:
                        st.markdown("### 🧪 SÉQUENCE TYPE")
                        st.markdown(f'<div class="dna-sequence">{gene_info["sequence_exemple"]}</div>', unsafe_allow_html=True)
                        
                        if st.button(f"📥 ANALYSER LA SÉQUENCE {gene_key}"):
                            st.session_state.current_sequence = gene_info["sequence_exemple"]
                            st.rerun()
                
                with gene_tabs[1]:
                    st.markdown("### 🧬 MUTATIONS CONNUES")
                    df_mutations = pd.DataFrame(gene_info['mutations'])
                    st.dataframe(df_mutations, use_container_width=True)
                    
                    # Graphique d'impact
                    impact_counts = Counter([m['impact'] for m in gene_info['mutations']])
                    fig = px.pie(values=list(impact_counts.values()), 
                               names=list(impact_counts.keys()),
                               title="Distribution des impacts des mutations")
                    st.plotly_chart(fig, use_container_width=True)
                
                with gene_tabs[2]:
                    st.markdown("### 📊 PHÉNOTYPES ASSOCIÉS")
                    for pheno in gene_info['phénotype']:
                        st.markdown(f"✅ {pheno}")
                    
                    st.markdown("### 🧬 QTL ASSOCIÉS")
                    for qtl in gene_info['qtl_associes']:
                        st.markdown(f"🔗 {qtl}")
                
                with gene_tabs[3]:
                    st.markdown("### 📚 PUBLICATIONS")
                    for pub in gene_info['publications']:
                        st.markdown(f"• {pub}")
                    
                    st.markdown("---")
                    st.markdown("### 💾 EXPORTER LES DONNÉES")
                    
                    col_export1, col_export2 = st.columns(2)
                    
                    with col_export1:
                        if st.button("📄 JSON"):
                            st.download_button(
                                label="Télécharger JSON",
                                data=json.dumps(gene_info, indent=2),
                                file_name=f"gene_{gene_key}_info.json",
                                mime="application/json"
                            )
                    
                    with col_export2:
                        if st.button("📊 CSV"):
                            # Convertir en format tabulaire
                            csv_data = pd.DataFrame([gene_info]).to_csv(index=False)
                            st.download_button(
                                label="Télécharger CSV",
                                data=csv_data,
                                file_name=f"gene_{gene_key}_info.csv",
                                mime="text/csv"
                            )
            else:
                st.warning(f"Gène '{gene_query}' non trouvé dans la base de données")
        
        # Liste complète des gènes
        with st.expander("📚 LISTE COMPLÈTE DES GÈNES OVINS (Cliquez pour développer)"):
            df_all_genes = pd.DataFrame([
                {
                    'Symbole': gene,
                    'Nom': info['nom_complet'],
                    'Chromosome': info['chromosome'],
                    'Fonction': info['fonction'][:100] + '...',
                    'Héritabilité': info['heritabilite']
                }
                for gene, info in genes_database.items()
            ])
            st.dataframe(df_all_genes, use_container_width=True)
    
    # Tab 3: GWAS & QTL
    with genetique_tabs[2]:
        st.markdown('<h2 class="section-header">📊 ANALYSE GWAS & QTL</h2>', unsafe_allow_html=True)
        
        # Interface GWAS
        col_gwas1, col_gwas2 = st.columns([2, 1])
        
        with col_gwas1:
            trait_gwas = st.selectbox(
                "**SÉLECTIONNEZ UN TRAIT POUR L'ANALYSE GWAS:**",
                ["Production laitière", "Teneur en matière grasse", "Croissance musculaire", 
                 "Fertilité", "Résistance aux maladies", "Taille de portée", "Longévité"]
            )
            
            significance_threshold = st.slider(
                "**SEUIL DE SIGNIFICATIVITÉ:**",
                min_value=1e-10,
                max_value=0.05,
                value=5e-8,
                format="%.2e"
            )
        
        with col_gwas2:
            st.markdown("""
            <div style='background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); 
                        color: white; padding: 20px; border-radius: 15px;'>
                <h4>⚙️ PARAMÈTRES GWAS</h4>
                <p>• Méthode: MLM</p>
                <p>• Covariables: Race, Âge</p>
                <p>• N SNPs: ~500,000</p>
                <p>• N individus: 1,000</p>
            </div>
            """, unsafe_allow_html=True)
        
        if st.button("🚀 LANCER L'ANALYSE GWAS", type="primary"):
            with st.spinner("📊 Exécution de l'analyse GWAS... Génération des Manhattan plots"):
                # Génération de données GWAS simulées
                np.random.seed(42)
                
                # Créer des données pour 27 chromosomes
                chromosomes = list(range(1, 28))
                all_data = []
                
                for chr_num in chromosomes:
                    n_snps = np.random.randint(50, 200)
                    positions = np.sort(np.random.choice(range(1, 1000000), n_snps, replace=False))
                    
                    # Générer des p-values avec quelques pics significatifs
                    p_values = np.random.exponential(0.1, n_snps)
                    
                    # Ajouter des pics artificiels pour certains chromosomes
                    if chr_num in [2, 6, 14]:  # Chromosomes avec QTL connus
                        peak_positions = np.random.choice(positions, 2)
                        for pos in peak_positions:
                            idx = np.where(positions == pos)[0][0]
                            p_values[idx] = np.random.uniform(1e-10, 1e-6)
                    
                    for pos, pval in zip(positions, p_values):
                        all_data.append({
                            'CHR': f'Chr{chr_num}',
                            'POS': pos,
                            'P': min(pval, 1.0)
                        })
                
                df_gwas = pd.DataFrame(all_data)
                df_gwas['-log10(P)'] = -np.log10(df_gwas['P'])
                
                # Manhattan plot
                st.markdown("### 📈 MANHATTAN PLOT")
                
                # Couleurs alternées pour les chromosomes
                colors = ['#2E7D32', '#4CAF50'] * 14
                chr_colors = {f'Chr{i+1}': colors[i % 2] for i in range(27)}
                
                fig = px.scatter(df_gwas, x='POS', y='-log10(P)', color='CHR',
                               color_discrete_map=chr_colors,
                               title=f"Manhattan Plot - GWAS pour {trait_gwas}",
                               labels={'POS': 'Position (bp)', '-log10(P)': '-log₁₀(p-value)'},
                               height=600)
                
                # Lignes de significativité
                fig.add_hline(y=-np.log10(significance_threshold), 
                            line_dash="dash", line_color="red",
                            annotation_text=f"Seuil: p < {significance_threshold:.2e}")
                
                fig.add_hline(y=-np.log10(0.05/len(df_gwas)), 
                            line_dash="dot", line_color="orange",
                            annotation_text="Bonferroni")
                
                st.plotly_chart(fig, use_container_width=True)
                
                # SNPs significatifs
                significant = df_gwas[df_gwas['P'] < significance_threshold]
                
                if len(significant) > 0:
                    st.success(f"✅ {len(significant)} SNPs significatifs détectés (p < {significance_threshold:.2e})")
                    
                    # Top 10 SNPs
                    top_snps = significant.nsmallest(10, 'P')
                    st.markdown("### 🏆 TOP 10 SNPs SIGNIFICATIFS")
                    st.dataframe(top_snps[['CHR', 'POS', 'P', '-log10(P)']])
                    
                    # Distribution par chromosome
                    chr_dist = significant['CHR'].value_counts().reset_index()
                    chr_dist.columns = ['Chromosome', 'Nombre de SNPs']
                    
                    fig2 = px.bar(chr_dist, x='Chromosome', y='Nombre de SNPs',
                                 title="Distribution des SNPs significatifs par chromosome")
                    st.plotly_chart(fig2, use_container_width=True)
                else:
                    st.info("ℹ️ Aucun SNP significatif détecté au seuil sélectionné")
                
                # QQ-plot
                st.markdown("### 📊 QQ-PLOT")
                
                # Générer des p-values attendues
                expected_p = np.sort(np.random.uniform(0, 1, len(df_gwas)))
                observed_p = np.sort(df_gwas['P'].values)
                
                fig3 = go.Figure()
                fig3.add_trace(go.Scatter(
                    x=-np.log10(expected_p),
                    y=-np.log10(observed_p),
                    mode='markers',
                    name='SNPs'
                ))
                
                # Ligne diagonale
                max_val = max(-np.log10(expected_p).max(), -np.log10(observed_p).max())
                fig3.add_trace(go.Scatter(
                    x=[0, max_val],
                    y=[0, max_val],
                    mode='lines',
                    name='Attendu',
                    line=dict(color='red', dash='dash')
                ))
                
                fig3.update_layout(
                    title="QQ-Plot - Distribution des p-values",
                    xaxis_title="-log₁₀(p-value attendu)",
                    yaxis_title="-log₁₀(p-value observé)",
                    showlegend=True
                )
                
                st.plotly_chart(fig3, use_container_width=True)
        
        # Section QTL
        st.markdown("---")
        st.markdown('<h2 class="section-header">🎯 BASE DE DONNÉES QTL OVINS</h2>', unsafe_allow_html=True)
        
        # Récupérer les QTL de la base de données
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM qtl_ovins")
        qtl_data = cursor.fetchall()
        
        if qtl_data:
            qtl_columns = ['ID', 'Nom', 'Chromosome', 'Start', 'End', 'Trait', 
                          'LOD', 'Variance', 'Peak', 'Genes', 'Race', 'Publication']
            
            df_qtl = pd.DataFrame(qtl_data, columns=qtl_columns)
            
            # Interface de filtrage
            col_filter1, col_filter2, col_filter3 = st.columns(3)
            
            with col_filter1:
                filter_trait = st.selectbox("Filtrer par trait:", 
                                          ["Tous"] + list(df_qtl['Trait'].unique()))
            
            with col_filter2:
                filter_chr = st.selectbox("Filtrer par chromosome:", 
                                        ["Tous"] + list(df_qtl['Chromosome'].unique()))
            
            with col_filter3:
                min_lod = st.slider("LOD minimum:", 0.0, 20.0, 5.0)
            
            # Appliquer les filtres
            filtered_df = df_qtl.copy()
            if filter_trait != "Tous":
                filtered_df = filtered_df[filtered_df['Trait'] == filter_trait]
            if filter_chr != "Tous":
                filtered_df = filtered_df[filtered_df['Chromosome'] == filter_chr]
            filtered_df = filtered_df[filtered_df['LOD'] >= min_lod]
            
            # Afficher les résultats
            st.dataframe(filtered_df[['Nom', 'Chromosome', 'Trait', 'LOD', 'Variance', 'Genes']], 
                        use_container_width=True)
            
            # Graphique des LOD scores
            fig_qtl = px.bar(filtered_df, x='Nom', y='LOD', color='Trait',
                           title="LOD Scores des QTL",
                           hover_data=['Variance', 'Genes'])
            st.plotly_chart(fig_qtl, use_container_width=True)
        else:
            st.info("Chargement des données QTL...")
    
    # Tab 4: Génétique des populations
    with genetique_tabs[3]:
        st.markdown('<h2 class="section-header">🧮 GÉNÉTIQUE DES POPULATIONS</h2>', unsafe_allow_html=True)
        
        pop_tabs = st.tabs(["Hardy-Weinberg", "Diversité", "F-statistiques", "Structure"])
        
        with pop_tabs[0]:
            st.markdown("### ⚖️ TEST D'ÉQUILIBRE DE HARDY-WEINBERG")
            
            # Interface de saisie
            genotypes_input = st.text_area(
                "**ENTREZ LES GÉNOTYPES** (un par ligne, format: AA, AB, BB, etc.):",
                height=150,
                placeholder="AA\nAB\nBB\nAA\nAB\nAA\nBB\nAB\nAA\nBB\nAB\nAA\nBB\nAA\nAB\nBB\nAA\nAB\nAA\nBB"
            )
            
            if st.button("📊 CALCULER HARDY-WEINBERG"):
                if genotypes_input:
                    genotypes = [g.strip().upper() for g in genotypes_input.split('\n') if g.strip()]
                    
                    hw_results = GeneticAnalyzerCustom.calculate_hardy_weinberg(genotypes)
                    
                    if hw_results:
                        # Résultats
                        col_hw1, col_hw2, col_hw3 = st.columns(3)
                        
                        with col_hw1:
                            st.markdown(f"""
                            <div class='metric-enhanced'>
                                <h3>χ²</h3>
                                <h2>{hw_results['chi_squared']}</h2>
                            </div>
                            """, unsafe_allow_html=True)
                        
                        with col_hw2:
                            st.markdown(f"""
                            <div class='metric-enhanced'>
                                <h3>p-value</h3>
                                <h2>{hw_results['p_value']:.6f}</h2>
                            </div>
                            """, unsafe_allow_html=True)
                        
                        with col_hw3:
                            status = "✅ ÉQUILIBRE" if hw_results['in_equilibrium'] else "⚠️ DÉSÉQUILIBRE"
                            st.markdown(f"""
                            <div class='metric-enhanced'>
                                <h3>STATUT</h3>
                                <h2>{status}</h2>
                            </div>
                            """, unsafe_allow_html=True)
                        
                        # Graphique observé vs attendu
                        genotypes_list = list(set(genotypes))
                        observed = [hw_results['observed'].get(g, 0) for g in genotypes_list]
                        expected = [hw_results['expected'].get(g, 0) for g in genotypes_list]
                        
                        fig_hw = go.Figure(data=[
                            go.Bar(name='Observé', x=genotypes_list, y=observed),
                            go.Bar(name='Attendu', x=genotypes_list, y=expected)
                        ])
                        
                        fig_hw.update_layout(
                            title="Distribution des génotypes - Observé vs Attendu",
                            barmode='group',
                            xaxis_title="Génotype",
                            yaxis_title="Fréquence"
                        )
                        
                        st.plotly_chart(fig_hw, use_container_width=True)
                        
                        # Fréquences alléliques
                        st.markdown("### 📊 FRÉQUENCES ALLÉLIQUES")
                        df_alleles = pd.DataFrame({
                            'Allèle': list(hw_results['allele_frequencies'].keys()),
                            'Fréquence': list(hw_results['allele_frequencies'].values())
                        })
                        
                        fig_alleles = px.bar(df_alleles, x='Allèle', y='Fréquence',
                                           title="Fréquences alléliques")
                        st.plotly_chart(fig_alleles, use_container_width=True)
        
        with pop_tabs[1]:
            st.markdown("### 🌍 ANALYSE DE DIVERSITÉ GÉNÉTIQUE")
            
            # Exemple de données
            st.info("Utilisez l'analyse Hardy-Weinberg ci-dessus pour générer des données de diversité")
            
            if 'hw_results' in locals() and hw_results:
                # Calculer la diversité
                allele_freq = hw_results['allele_frequencies']
                heterozygosity_obs = sum(1 for g in genotypes if g[0] != g[1]) / len(genotypes)
                heterozygosity_exp = 1 - sum(f**2 for f in allele_freq.values())
                
                # Mesures de diversité
                col_div1, col_div2, col_div3 = st.columns(3)
                
                with col_div1:
                    st.metric("Ho (Hétérozygotie observée)", f"{heterozygosity_obs:.4f}")
                
                with col_div2:
                    st.metric("He (Hétérozygotie attendue)", f"{heterozygosity_exp:.4f}")
                
                with col_div3:
                    fis = 1 - (heterozygosity_obs / heterozygosity_exp) if heterozygosity_exp > 0 else 0
                    st.metric("Fis (Consanguinité)", f"{fis:.4f}")
                
                # Indice de Shannon
                shannon = -sum(f * math.log(f) for f in allele_freq.values() if f > 0)
                st.metric("Indice de Shannon", f"{shannon:.4f}")
    
    # Tab 5: Pedigrees
    with genetique_tabs[4]:
        st.markdown('<h2 class="section-header">🌳 ANALYSE DE PEDIGREES</h2>', unsafe_allow_html=True)
        
        pedigree_input = st.text_area(
            "**ENTREZ LE PEDIGREE** (Format: Animal,Père,Mère):",
            height=200,
            placeholder="Animal1,,\\nAnimal2,Père1,Mère1\\nAnimal3,Père2,Mère2\\nAnimal4,Père1,Mère2\\nAnimal5,Père3,Mère3"
        )
        
        if st.button("🌳 ANALYSER LE PEDIGREE"):
            if pedigree_input:
                # Parser le pedigree
                pedigree_data = []
                for line in pedigree_input.strip().split('\n'):
                    parts = [p.strip() for p in line.split(',')]
                    if len(parts) >= 3:
                        pedigree_data.append((parts[0], parts[1] if parts[1] else None, 
                                            parts[2] if parts[2] else None))
                
                if pedigree_data:
                    # Calculs basiques
                    animals = [p[0] for p in pedigree_data]
                    sires = [p[1] for p in pedigree_data if p[1]]
                    dams = [p[2] for p in pedigree_data if p[2]]
                    
                    col_ped1, col_ped2, col_ped3 = st.columns(3)
                    
                    with col_ped1:
                        st.metric("Animaux", len(animals))
                    
                    with col_ped2:
                        st.metric("Pères uniques", len(set(sires)))
                    
                    with col_ped3:
                        st.metric("Mères uniques", len(set(dams)))
                    
                    # Statistiques de parenté
                    st.markdown("### 📊 STATISTIQUES DE PARENTÉ")
                    
                    # Calcul simplifié du coefficient de consanguinité
                    inbreeding_data = []
                    for animal, sire, dam in pedigree_data:
                        if sire and dam:
                            # Pour la démo, coefficient aléatoire
                            coeff = round(random.uniform(0.0, 0.25), 4)
                            inbreeding_data.append({
                                'Animal': animal,
                                'Père': sire,
                                'Mère': dam,
                                'Coefficient': coeff
                            })
                    
                    if inbreeding_data:
                        df_inbreeding = pd.DataFrame(inbreeding_data)
                        st.dataframe(df_inbreeding.sort_values('Coefficient', ascending=False))
                        
                        # Graphique
                        fig_ped = px.histogram(df_inbreeding, x='Coefficient',
                                              title="Distribution des coefficients de consanguinité",
                                              nbins=20)
                        st.plotly_chart(fig_ped, use_container_width=True)
    
    # Tab 6: Import/Export
    with genetique_tabs[5]:
        st.markdown('<h2 class="section-header">💾 IMPORTATION & EXPORTATION DE DONNÉES</h2>', unsafe_allow_html=True)
        
        format_tabs = st.tabs(["VCF", "PLINK", "FASTA", "EXCEL"])
        
        with format_tabs[0]:
            st.markdown("### 📄 FORMAT VCF (Variant Call Format)")
            
            # Générer un VCF exemple
            vcf_example = """##fileformat=VCFv4.2
##fileDate=20240115
##source=OvinManagerPro_v3.0
##reference=Ovis_aries_4.0
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOuled_Djellal\tRazè\tHamra
2\t123456\trs123456\tA\tG\t100\tPASS\tAF=0.42;DP=50\tGT\t0/0\t0/1\t1/1
6\t654321\trs654321\tC\tT\t150\tPASS\tAF=0.18;DP=60\tGT\t0/1\t0/0\t1/1
14\t789123\trs789123\tG\tA\t200\tPASS\tAF=0.33;DP=45\tGT\t1/1\t0/1\t0/0
5\t321654\trs321654\tT\tC\t180\tPASS\tAF=0.25;DP=55\tGT\t0/0\t1/1\t0/1"""
            
            st.code(vcf_example, language='vcf')
            
            st.download_button(
                label="📥 TÉLÉCHARGER EXEMPLE VCF",
                data=vcf_example,
                file_name="ovin_genotypes_example.vcf",
                mime="text/plain"
            )
        
        with format_tabs[1]:
            st.markdown("### 📊 FORMAT PLINK (.bed/.bim/.fam)")
            st.info("Génération de fichiers PLINK en cours de développement...")
            
            # Exemple de .fam
            fam_example = """Ouled_Djellal_001  Ouled_Djellal_001  0 0 1 1
Razè_001         Razè_001         0 0 2 1
Hamra_001        Hamra_001        0 0 1 1
Dman_001         Dman_001         0 0 2 1
Saharienne_001   Saharienne_001   0 0 1 1"""
            
            st.code(fam_example, language='text')
        
        with format_tabs[2]:
            st.markdown("### 🧬 FORMAT FASTA")
            
            fasta_input = st.text_area("Entrez des séquences au format FASTA:", height=200)
            
            if fasta_input:
                st.code(fasta_input, language='fasta')
        
        with format_tabs[3]:
            st.markdown("### 📈 EXPORT EXCEL")
            
            # Exporter les données de la base
            cursor = conn.cursor()
            
            tables = ['ovins', 'genetic_markers', 'qtl_ovins']
            
            for table in tables:
                cursor.execute(f"SELECT * FROM {table} LIMIT 10")
                data = cursor.fetchall()
                
                if data:
                    st.markdown(f"#### Table: {table}")
                    
                    # Récupérer les noms des colonnes
                    cursor.execute(f"PRAGMA table_info({table})")
                    columns = [col[1] for col in cursor.fetchall()]
                    
                    df_table = pd.DataFrame(data, columns=columns)
                    st.dataframe(df_table)
                    
                    # Bouton d'export
                    excel_data = df_table.to_csv(index=False)
                    st.download_button(
                        label=f"📥 Exporter {table} (CSV)",
                        data=excel_data,
                        file_name=f"{table}_export.csv",
                        mime="text/csv"
                    )

# ========== NAVIGATION PRINCIPALE ==========

# Sidebar
with st.sidebar:
    st.markdown("""
    <div style='text-align: center; padding: 20px; background: linear-gradient(135deg, #2E7D32 0%, #4CAF50 100%); 
                color: white; border-radius: 10px; margin-bottom: 20px;'>
        <h2>🧬 OVIN MANAGER PRO</h2>
        <p>Version Génétique Avancée</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("### 📍 NAVIGATION GÉNÉTIQUE")
    
    page = st.radio(
        "LABORATOIRE DE GÉNÉTIQUE",
        ["🧬 GÉNÉTIQUE AVANCÉE", 
         "📊 GESTION TROUPEAU", 
         "🥛 ANALYSE LAIT", 
         "🤰 GESTATION",
         "📈 STATISTIQUES",
         "⚙️ PARAMÈTRES"]
    )
    
    st.markdown("---")
    
    st.markdown("### 📊 STATISTIQUES GÉNÉTIQUES")
    
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM genetic_markers")
    markers = cursor.fetchone()[0]
    st.metric("Marqueurs", f"{markers:,}")
    
    cursor.execute("SELECT COUNT(*) FROM qtl_ovins")
    qtls = cursor.fetchone()[0]
    st.metric("QTL", qtls)
    
    cursor.execute("SELECT COUNT(DISTINCT gene_symbol) FROM genetic_markers")
    genes = cursor.fetchone()[0]
    st.metric("Gènes", genes)

# Affichage de la page principale
if page == "🧬 GÉNÉTIQUE AVANCÉE":
    page_genetique_pro()
elif page == "📊 GESTION TROUPEAU":
    st.info("Module de gestion du troupeau - À intégrer")
elif page == "🥛 ANALYSE LAIT":
    st.info("Module d'analyse laitière - À intégrer")
elif page == "🤰 GESTION":
    st.info("Module de gestation - À intégrer")
elif page == "📈 STATISTIQUES":
    st.info("Module de statistiques - À intégrer")
elif page == "⚙️ PARAMÈTRES":
    st.info("Module de paramètres - À intégrer")

# Pied de page
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🧬 <strong>OVIN MANAGER PRO - LABORATOIRE DE GÉNÉTIQUE AVANCÉE</strong></p>
    <p>Version 3.0 | Pour la recherche et la sélection génétique ovine</p>
    <p>© 2024 - Tous droits réservés</p>
</div>
""", unsafe_allow_html=True)
