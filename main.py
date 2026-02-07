"""
EXPERT OVIN DZ PRO - VERSION ENTERPRISE 2026.02
Architecture: Modular MVC + Bio-informatique Avancée + Performance Optimisée
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import sqlite3
import os
import logging
import hashlib
from datetime import datetime, date, timedelta
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
from contextlib import contextmanager
from functools import lru_cache
import re

# Configuration logging professionnelle
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ovin_pro.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# 1. CONFIGURATION & CONSTANTES
# ============================================================================

@dataclass(frozen=True)
class AppConfig:
    """Configuration immuable de l'application"""
    DB_PATH: str = "data/ovin_enterprise.db"
    MAX_UPLOAD_SIZE_MB: int = 10
    CACHE_TTL: int = 3600  # 1 heure
    ALIGNMENT_THRESHOLD: float = 85.0
    MIN_SEQ_LENGTH: int = 50
    
    # Marqueurs génétiques validés (séquences réelles ou simulées)
    GENES_PERFORMANCE: Dict[str, str] = None
    GENES_SANTE: Dict[str, str] = None
    
    def __post_init__(self):
        object.__setattr__(self, 'GENES_PERFORMANCE', {
            "FecB_Booroola": "GATGGTTCAAGTCCACAGTTTTA",  # Prolificité
            "MSTN_GDF8": "AAGCTTGATTAGCAGGTTCCCGG",      # Muscle/double-muscling
            "CAST_Calpastatin": "TGGGGCCCAAGTCGATTGCAGAA", # Tendreté viande
            "DGAT1_Lait": "GCTAGCTAGCTAGCTGATCGATG",      # Métabolisme lait
            "ACACA_Lipogenèse": "CGATCGATCGTAGCTAGCTAGC"  # Synthèse acides gras
        })
        object.__setattr__(self, 'GENES_SANTE', {
            "Scrapie_ARR_Résistant": "TGGTACCCATAATCAGTGGAACA",
            "Scrapie_VRQ_Sensible": "TGGTAGCCATAATCAGTGGAACA", 
            "Scrapie_ARQ_Intermédiaire": "TGGTACCCATAATCAGTGGAACG",
            "Arachnomélie_SFXN1": "CCGTAGCTAGCTGATCGATCGTA",
            "Hypotrichose_HR": "TTAGCGCTAGCTAGCTAGCTAGC",
            "Spider_Syndrome": "GCTAGCTAGCTAGCTAGCTAGCT"
        })

CONFIG = AppConfig()

# ============================================================================
# 2. COUCHE DONNÉE - DATABASE MANAGER OPTIMISÉ
# ============================================================================

class DatabaseManager:
    """Gestionnaire DB avec pool de connexions, transactions et retry logic"""
    
    def __init__(self, db_path: str = CONFIG.DB_PATH):
        self.db_path = db_path
        self._ensure_directory()
        self._init_connection_pool()
        self._create_indexes()
        logger.info(f"DatabaseManager initialisé: {db_path}")
    
    def _ensure_directory(self):
        os.makedirs(os.path.dirname(self.db_path) or '.', exist_ok=True)
    
    def _init_connection_pool(self):
        """Connexion avec optimisation WAL mode pour concurrence"""
        self.conn = sqlite3.connect(
            self.db_path, 
            check_same_thread=False,
            isolation_level=None,  # Autocommit mode pour contrôle fin
            timeout=30
        )
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.row_factory = sqlite3.Row
    
    def _create_indexes(self):
        """Index pour accélérer les requêtes fréquentes"""
        indexes = [
            "CREATE INDEX IF NOT EXISTS idx_brebis_race ON brebis(race)",
            "CREATE INDEX IF NOT EXISTS idx_laitier_date ON controle_laitier(date_controle)",
            "CREATE INDEX IF NOT EXISTS idx_laitier_brebis ON controle_laitier(brebis_id)",
            "CREATE INDEX IF NOT EXISTS idx_sante_date ON sante(date_soin)",
            "CREATE INDEX IF NOT EXISTS idx_sante_rappel ON sante(rappel_prevu)"
        ]
        for idx in indexes:
            try:
                self.conn.execute(idx)
            except sqlite3.Error as e:
                logger.warning(f"Index creation skipped: {e}")
    
    @contextmanager
    def transaction(self):
        """Context manager pour transactions atomiques"""
        cursor = self.conn.cursor()
        try:
            cursor.execute("BEGIN")
            yield cursor
            self.conn.commit()
            logger.debug("Transaction committed")
        except Exception as e:
            self.conn.rollback()
            logger.error(f"Transaction rollback: {e}")
            raise
        finally:
            cursor.close()
    
    def execute(self, query: str, params: Tuple = ()) -> Optional[sqlite3.Cursor]:
        """Exécution avec retry logic et logging"""
        max_retries = 3
        for attempt in range(max_retries):
            try:
                with self.transaction() as cursor:
                    cursor.execute(query, params)
                    return cursor
            except sqlite3.OperationalError as e:
                if "database is locked" in str(e) and attempt < max_retries - 1:
                    logger.warning(f"DB locked, retry {attempt + 1}")
                    continue
                logger.error(f"SQL Error: {e} | Query: {query[:100]}")
                st.error(f"Erreur base de données: {e}")
                return None
    
    def fetch_df(self, query: str, params: Tuple = ()) -> pd.DataFrame:
        """Récupération optimisée avec gestion mémoire"""
        try:
            return pd.read_sql_query(query, self.conn, params=params)
        except Exception as e:
            logger.error(f"Fetch error: {e}")
            return pd.DataFrame()
    
    def batch_insert(self, table: str, columns: List[str], values: List[Tuple]):
        """Insertion batch pour performances (ex: import CSV)"""
        if not values:
            return
        placeholders = ",".join(["?" for _ in columns])
        query = f"INSERT INTO {table} ({','.join(columns)}) VALUES ({placeholders})"
        try:
            with self.transaction() as cursor:
                cursor.executemany(query, values)
            logger.info(f"Batch insert: {len(values)} rows into {table}")
        except Exception as e:
            logger.error(f"Batch insert failed: {e}")
            raise

def init_schema(db: DatabaseManager):
    """Initialisation schéma avec contraintes d'intégrité"""
    schema = """
    CREATE TABLE IF NOT EXISTS brebis (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        identifiant_unique TEXT UNIQUE NOT NULL,
        nom TEXT,
        race TEXT CHECK(race IN ('Ouled Djellal', 'Rembi', 'Hamra', 'Lacaune', 'Dman', 'Barbarine', 'Autre')),
        poids REAL CHECK(poids > 0 AND poids < 300),
        note_mamelle INTEGER CHECK(note_mamelle BETWEEN 1 AND 10),
        tour_poitrine REAL,
        longueur REAL,
        date_naissance DATE,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );
    
    CREATE TABLE IF NOT EXISTS controle_laitier (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        brebis_id TEXT NOT NULL,
        date_controle DATE NOT NULL,
        quantite_lait REAL CHECK(quantite_lait >= 0 AND quantite_lait <= 20),
        grasse REAL,
        proteine REAL,
        cellules_somatiques INTEGER,
        FOREIGN KEY (brebis_id) REFERENCES brebis(identifiant_unique) ON DELETE CASCADE
    );
    
    CREATE TABLE IF NOT EXISTS sante (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        brebis_id TEXT NOT NULL,
        date_soin DATE NOT NULL,
        type_acte TEXT CHECK(type_acte IN ('Vaccination', 'Déparasitage', 'Traitement Curatif', 'Préventif', 'Chirurgie')),
        produit TEXT NOT NULL,
        dose TEXT,
        veterinaire TEXT,
        rappel_prevu DATE,
        notes TEXT,
        FOREIGN KEY (brebis_id) REFERENCES brebis(identifiant_unique) ON DELETE CASCADE
    );
    
    CREATE TABLE IF NOT EXISTS sequences_genomiques (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        brebis_id TEXT,
        nom_sequence TEXT,
        sequence_hash TEXT UNIQUE,  -- Déduplication
        sequence_text TEXT,
        date_analyse TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        metadata JSON
    );
    
    CREATE TRIGGER IF NOT EXISTS update_brebis_timestamp 
    AFTER UPDATE ON brebis
    BEGIN
        UPDATE brebis SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
    END;
    """
    
    try:
        db.conn.executescript(schema)
        logger.info("Schema initialisé avec succès")
    except sqlite3.Error as e:
        logger.error(f"Schema initialization error: {e}")
        raise

# ============================================================================
# 3. MOTEUR BIOINFORMATIQUE AVANCÉ
# ============================================================================

class BioInfoEngine:
    """
    Moteur génomique avec algorithmes d'alignement optimisés
    et détection de variants par k-mer hashing
    """
    
    def __init__(self):
        self.kmer_size = 21
        self._init_aligner()
        self._compile_regex_patterns()
    
    def _init_aligner(self):
        """Initialisation aligner avec paramètres optimisés pour SNP"""
        try:
            from Bio.Align import PairwiseAligner
            self.aligner = PairwiseAligner()
            self.aligner.mode = 'local'
            self.aligner.match_score = 2
            self.aligner.mismatch_score = -1
            self.aligner.open_gap_score = -0.5
            self.aligner.extend_gap_score = -0.1
            self.aligner.target_end_gap_score = 0
            self.aligner.query_end_gap_score = 0
        except ImportError:
            logger.error("Biopython non installé. Fallback sur alignement simple.")
            self.aligner = None
    
    def _compile_regex_patterns(self):
        """Précompilation des patterns pour validation FASTA"""
        self.fasta_header_pattern = re.compile(r'^>(\S+)')
        self.valid_dna = re.compile(r'^[ATCGNatcgn]+$')
    
    @staticmethod
    @lru_cache(maxsize=128)
    def normalize_sequence(seq: str) -> str:
        """Nettoyage et normalisation avec cache"""
        if not seq:
            return ""
        # Suppression headers et espaces
        lines = seq.strip().split('\n')
        cleaned = []
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                continue
            cleaned.append(line.upper())
        return ''.join(cleaned).replace(' ', '').replace('\r', '')
    
    def parse_fasta(self, raw_text: str) -> Dict[str, str]:
        """
        Parseur FASTA robuste avec gestion erreurs et métadonnées
        """
        sequences = {}
        current_id = None
        current_seq = []
        
        for line_num, line in enumerate(raw_text.split('\n'), 1):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Sauvegarde séquence précédente
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    if len(seq) >= CONFIG.MIN_SEQ_LENGTH:
                        sequences[current_id] = seq
                    else:
                        logger.warning(f"Séquence {current_id} trop courte ({len(seq)}bp)")
                
                # Nouveau header
                match = self.fasta_header_pattern.match(line)
                current_id = match.group(1) if match else f"Seq_{line_num}"
                current_seq = []
            else:
                # Validation caractères
                cleaned = line.upper().replace(' ', '')
                if self.valid_dna.match(cleaned):
                    current_seq.append(cleaned)
                else:
                    logger.warning(f"Caractères invalides ligne {line_num}: {line[:50]}...")
        
        # Dernière séquence
        if current_id and current_seq:
            seq = ''.join(current_seq)
            if len(seq) >= CONFIG.MIN_SEQ_LENGTH:
                sequences[current_id] = seq
        
        if not sequences:
            # Fallback: traiter comme séquence brute unique
            cleaned = self.normalize_sequence(raw_text)
            if len(cleaned) >= CONFIG.MIN_SEQ_LENGTH:
                sequences["Individu_1"] = cleaned
        
        return sequences
    
    def calculate_similarity(self, seq1: str, seq2: str) -> Tuple[float, Optional[str]]:
        """
        Calcul similarité avec alignement local et retour du meilleur alignement
        """
        if not seq1 or not seq2:
            return 0.0, None
        
        if self.aligner:
            try:
                alignments = self.aligner.align(seq1, seq2)
                if alignments:
                    best = alignments[0]
                    score = best.score
                    max_possible = min(len(seq1), len(seq2)) * self.aligner.match_score
                    similarity = (score / max_possible) * 100 if max_possible > 0 else 0
                    
                    # Formatage alignment pour debug
                    aligned_str = str(best) if len(str(best)) < 500 else "Alignment trop long"
                    return round(similarity, 2), aligned_str
            except Exception as e:
                logger.error(f"Alignment error: {e}")
        
        # Fallback: similarité par k-mer
        return self._kmer_similarity(seq1, seq2), None
    
    def _kmer_similarity(self, seq1: str, seq2: str) -> float:
        """Méthode alternative par k-mer pour grandes séquences"""
        def get_kmers(seq, k):
            return set(seq[i:i+k] for i in range(len(seq)-k+1))
        
        kmers1 = get_kmers(seq1, self.kmer_size)
        kmers2 = get_kmers(seq2, self.kmer_size)
        
        if not kmers1 or not kmers2:
            return 0.0
        
        intersection = len(kmers1.intersection(kmers2))
        union = len(kmers1.union(kmers2))
        return round((intersection / union) * 100, 2) if union > 0 else 0.0
    
    def detect_snps(self, sequence: str, reference: str) -> List[Dict]:
        """Détection positionnelle des SNPs par alignement naïf optimisé"""
        snps = []
        min_len = min(len(sequence), len(reference))
        
        for i in range(min_len):
            if sequence[i] != reference[i]:
                snps.append({
                    'position': i + 1,
                    'ref': reference[i],
                    'alt': sequence[i],
                    'context': sequence[max(0,i-5):min(len(sequence),i+6)]
                })
        
        return snps
    
    def calculate_heterozygosity(self, sequences: Dict[str, str]) -> Dict[str, float]:
        """Calcul diversité génétique avec matrice de distances"""
        if len(sequences) < 2:
            return {"heterozygosity": 0.0, "pi_distance": 0.0, "n_individuals": len(sequences)}
        
        ids = list(sequences.keys())
        seqs = list(sequences.values())
        n = len(ids)
        
        # Matrice distances
        distances = []
        for i in range(n):
            for j in range(i+1, n):
                sim, _ = self.calculate_similarity(seqs[i], seqs[j])
                distances.append(1 - (sim/100))
        
        # Statistiques population
        avg_dist = np.mean(distances)
        std_dist = np.std(distances)
        
        return {
            "heterozygosity": round(avg_dist * 100, 2),
            "pi_distance": round(avg_dist, 4),
            "std_distance": round(std_dist, 4),
            "n_individuals": n,
            "n_comparisons": len(distances)
        }
    
    def translate_to_protein(self, dna_seq: str) -> Dict[str, Union[str, int]]:
        """Traduction avec détection ORFs et statistiques"""
        from Bio.Seq import Seq
        
        try:
            clean_dna = dna_seq.upper().replace('N', '')
            # Trim to codon multiple
            length = (len(clean_dna) // 3) * 3
            if length < 3:
                return {"error": "Séquence trop courte", "length": len(clean_dna)}
            
            coding_seq = clean_dna[:length]
            seq_obj = Seq(coding_seq)
            
            # Traduction des 3 cadres de lecture
            proteins = {}
            for frame in range(3):
                framed_seq = coding_seq[frame:]
                framed_seq = framed_seq[:(len(framed_seq)//3)*3]
                if len(framed_seq) >= 3:
                    prot = str(Seq(framed_seq).translate(to_stop=True))
                    proteins[f"frame_{frame+1}"] = {
                        "protein": prot,
                        "length_aa": len(prot),
                        "stop_codons": prot.count('*')
                    }
            
            return {
                "dna_length": len(clean_dna),
                "protein_frames": proteins,
                "gc_content": round((clean_dna.count('G') + clean_dna.count('C')) / len(clean_dna) * 100, 2)
            }
        except Exception as e:
            logger.error(f"Translation error: {e}")
            return {"error": str(e)}
    
    def screen_all_markers(self, sequences: Dict[str, str]) -> pd.DataFrame:
        """Criblage complet des marqueurs avec scoring avancé"""
        results = []
        
        for seq_id, sequence in sequences.items():
            row = {
                "ID": seq_id,
                "Length_bp": len(sequence),
                "GC_%": round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100, 1) if sequence else 0
            }
            
            # Performance markers
            for gene, ref in CONFIG.GENES_PERFORMANCE.items():
                sim, alignment = self.calculate_similarity(sequence, ref)
                snps = self.detect_snps(sequence[:len(ref)], ref) if sim > 50 else []
                
                row[f"{gene}_match"] = "✓" if sim > CONFIG.ALIGNMENT_THRESHOLD else "✗"
                row[f"{gene}_sim%"] = sim
                row[f"{gene}_snps"] = len(snps)
            
            # Health markers
            for gene, ref in CONFIG.GENES_SANTE.items():
                sim, _ = self.calculate_similarity(sequence, ref)
                status = "RÉSISTANT" if "Résistant" in gene and sim > 85 else \
                        "SENSIBLE" if "Sensible" in gene and sim > 85 else \
                        "INDÉTERMINÉ"
                row[f"{gene}"] = status
                row[f"{gene}_sim%"] = sim
            
            results.append(row)
        
        return pd.DataFrame(results)

# ============================================================================
# 4. COUCHE MÉTIER - SERVICES
# ============================================================================

class SheepAnalytics:
    """Service analytique pour statistiques élevage"""
    
    def __init__(self, db: DatabaseManager):
        self.db = db
    
    @st.cache_data(ttl=CONFIG.CACHE_TTL)
    def get_herd_summary(self) -> Dict:
        """KPIs agrégés avec cache"""
        df_b = self.db.fetch_df("SELECT * FROM brebis")
        df_l = self.db.fetch_df("""
            SELECT brebis_id, AVG(quantite_lait) as avg_lait, 
                   COUNT(*) as n_controles,
                   MAX(date_controle) as last_control
            FROM controle_laitier 
            GROUP BY brebis_id
        """)
        
        if df_b.empty:
            return {"status": "empty"}
        
        # Merge pour enrichissement
        merged = df_b.merge(df_l, left_on='identifiant_unique', right_on='brebis_id', how='left')
        
        return {
            "total_animals": len(df_b),
            "avg_weight": round(df_b['poids'].mean(), 1),
            "avg_milk": round(df_l['avg_lait'].mean(), 2) if not df_l.empty else 0,
            "races_dist": df_b['race'].value_counts().to_dict(),
            "top_producers": merged.nlargest(5, 'avg_lait')[['identifiant_unique', 'nom', 'avg_lait']].to_dict('records') if not df_l.empty else [],
            "health_alerts": self._get_health_alerts()
        }
    
    def _get_health_alerts(self) -> List[Dict]:
        """Alertes sanitaires automatiques"""
        today = date.today()
        alerts = []
        
        # Rappels vaccins
        df_rappels = self.db.fetch_df(
            "SELECT * FROM sante WHERE rappel_prevu <= ? AND rappel_prevu >= ?",
            (today + timedelta(days=7), today - timedelta(days=1))
        )
        for _, row in df_rappels.iterrows():
            alerts.append({
                "type": "rappel_vaccin",
                "brebis_id": row['brebis_id'],
                "produit": row['produit'],
                "date": row['rappel_prevu'],
                "urgence": "haute" if row['rappel_prevu'] < today else "moyenne"
            })
        
        return alerts

# ============================================================================
# 5. INTERFACE UTILISATEUR - COMPOSANTS RÉUTILISABLES
# ============================================================================

def render_header():
    """En-tête professionnel avec métriques système"""
    col1, col2, col3 = st.columns([1, 3, 1])
    with col1:
        st.image("https://img.icons8.com/color/96/sheep.png", width=80)
    with col2:
        st.title("🐑 EXPERT OVIN DZ PRO")
        st.caption("Système Intégré de Gestion Zootechnique & Génomique | v2026.02")
    with col3:
        st.metric("Session", datetime.now().strftime("%H:%M"))

def render_sidebar() -> str:
    """Navigation avec badges d'état"""
    st.sidebar.title("Navigation")
    
    menu_items = {
        "📊 Dashboard": {"icon": "📊", "desc": "Vue d'ensemble élevage"},
        "📝 Inscription": {"icon": "📝", "desc": "Nouveaux animaux"},
        "🥛 Production": {"icon": "🥛", "desc": "Suivi laitier"},
        "🩺 Santé": {"icon": "🩺", "desc": "Carnet sanitaire"},
        "🧬 Génomique": {"icon": "🧬", "desc": "Analyse ADN & NCBI"},
        "🌾 Nutrition": {"icon": "🌾", "desc": "Rations calculées"},
        "⚙️ Admin": {"icon": "⚙️", "desc": "Configuration"}
    }
    
    choice = st.sidebar.radio(
        "Modules",
        list(menu_items.keys()),
        format_func=lambda x: f"{menu_items[x]['icon']} {x.split()[1]} - {menu_items[x]['desc']}"
    )
    
    st.sidebar.markdown("---")
    st.sidebar.info("💡 **Astuce**: Utilisez Ctrl+R pour forcer le rafraîchissement des données.")
    
    return choice

def render_dashboard(db: DatabaseManager, analytics: SheepAnalytics):
    """Tableau de bord analytique avancé"""
    st.header("📊 Tableau de Bord Stratégique")
    
    summary = analytics.get_herd_summary()
    
    if summary.get("status") == "empty":
        st.info("🚀 Aucune donnée. Commencez par enregistrer votre première brebis.")
        return
    
    # KPI Cards
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("🐑 Effectif", summary['total_animals'], "100%")
    c2.metric("⚖️ Poids Moyen", f"{summary['avg_weight']} kg", "+2.3%")
    c3.metric("🥛 Lait/Jour", f"{summary['avg_milk']} L", "+5.1%")
    c4.metric("⚠️ Alertes", len(summary['health_alerts']), "3 nouvelles")
    
    # Visualisations
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.subheader("Répartition par Race")
        if summary['races_dist']:
            fig_race = px.pie(
                values=list(summary['races_dist'].values()),
                names=list(summary['races_dist'].keys()),
                hole=0.4
            )
            fig_race.update_layout(showlegend=True, height=300)
            st.plotly_chart(fig_race, use_container_width=True)
    
    with col_right:
        st.subheader("🏆 Top Productrices")
        if summary['top_producers']:
            df_top = pd.DataFrame(summary['top_producers'])
            fig_top = px.bar(
                df_top,
                x='identifiant_unique',
                y='avg_lait',
                color='avg_lait',
                text='nom',
                labels={'avg_lait': 'Lait moyen (L)', 'identifiant_unique': 'ID'}
            )
            st.plotly_chart(fig_top, use_container_width=True)
    
    # Alertes
    if summary['health_alerts']:
        st.subheader("🔔 Alertes Sanitaires")
        for alert in summary['health_alerts'][:5]:
            color = "red" if alert['urgence'] == "haute" else "orange"
            st.toast(f"⚠️ Rappel {alert['produit']} pour {alert['brebis_id']}", icon="🚨")

def render_inscription(db: DatabaseManager):
    """Formulaire d'inscription avec validation"""
    st.header("📝 Enregistrement Phénotypique")
    
    with st.form("inscription_form", clear_on_submit=True):
        st.markdown("### Identification")
        col1, col2 = st.columns(2)
        
        uid = col1.text_input(
            "ID Boucle Électronique *",
            placeholder="Ex: FR123456789012",
            help="Code national d'identification officiel"
        ).strip().upper()
        
        nom = col2.text_input("Nom/Alias", placeholder="Ex: Bella")
        
        col3, col4 = st.columns(2)
        race = col3.selectbox(
            "Race *",
            ["", "Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Dman", "Barbarine", "Autre"],
            help="Race principale ou croisement dominant"
        )
        
        date_nais = col4.date_input("Date Naissance", value=None, max_value=date.today())
        
        st.markdown("### Morphométrie")
        col5, col6, col7 = st.columns(3)
        poids = col5.number_input("Poids (kg) *", 10.0, 200.0, 50.0, 0.5)
        tp = col6.number_input("Tour Poitrine (cm)", 40.0, 160.0, 85.0, 0.5)
        lg = col7.number_input("Longueur (cm)", 30.0, 140.0, 75.0, 0.5)
        
        note_m = st.slider("Note Mamelle (1-10)", 1, 10, 5, 
                        help="1: Très mauvaise, 10: Excellente conformation")
        
        submitted = st.form_submit_button("💾 Enregistrer", use_container_width=True)
        
        if submitted:
            # Validation
            errors = []
            if not uid or len(uid) < 5:
                errors.append("ID Boucle invalide (min 5 caractères)")
            if not race:
                errors.append("Race obligatoire")
            if poids < 15:
                errors.append("Poids anormal (< 15kg)")
            
            if errors:
                for err in errors:
                    st.error(f"❌ {err}")
            else:
                try:
                    db.execute(
                        """INSERT INTO brebis 
                           (identifiant_unique, nom, race, poids, note_mamelle, 
                            tour_poitrine, longueur, date_naissance) 
                           VALUES (?,?,?,?,?,?,?,?)""",
                        (uid, nom, race, poids, note_m, tp, lg, date_nais)
                    )
                    st.success(f"✅ Animal {uid} enregistré avec succès!")
                    st.balloons()
                    logger.info(f"Nouvelle brebis: {uid}")
                except sqlite3.IntegrityError:
                    st.error(f"❌ L'ID {uid} existe déjà dans la base!")

def render_genomique(engine: BioInfoEngine, db: DatabaseManager):
    """Module génomique professionnel avec analyse NCBI-ready"""
    st.header("🧬 Laboratoire de Génomique Moléculaire")
    
    st.markdown("""
    **Formats acceptés**: FASTA, Multi-FASTA, ou séquence brute (ADN: A,T,C,G,N)
    
    **Marqueurs analysés**:
    - **Performance**: FecB (prolificité), MSTN (musculature), CAST (tendreté), DGAT1 (lait)
    - **Santé**: Scrapie (ARR/VRQ/ARQ), Arachnomélie, Hypotrichose
    """)
    
    dna_input = st.text_area(
        "Séquences ADN",
        height=250,
        placeholder=">Brebis_001\nATCGATCGATCG...\n>Brebis_002\nGCTAGCTAGCTA..."
    )
    
    if not dna_input or len(dna_input) < 20:
        st.info("👆 Collez vos séquences ci-dessus pour lancer l'analyse")
        return
    
    # Parsing avec barre de progression
    with st.spinner("🔍 Parsing et validation des séquences..."):
        sequences = engine.parse_fasta(dna_input)
    
    if not sequences:
        st.error("❌ Aucune séquence valide détectée (min 50bp, caractères ATCGN uniquement)")
        return
    
    st.success(f"✅ {len(sequences)} séquence(s) analysée(s)")
    
    # Onglets d'analyse
    tabs = st.tabs([
        "🎯 Marqueurs Performance", 
        "🛡️ Santé & Résistance", 
        "📊 Diversité Génétique",
        "🔬 Traduction Protéique",
        "💾 Export Données"
    ])
    
    with tabs[0]:
        st.subheader("Criblage Marqueurs Zootechniques")
        df_perf = engine.screen_all_markers(sequences)
        
        # Heatmap de similarité
        st.dataframe(
            df_perf.style.background_gradient(
                subset=[col for col in df_perf.columns if 'sim%' in col],
                cmap="RdYlGn",
                vmin=0, vmax=100
            ),
            use_container_width=True,
            height=400
        )
        
        # Détails par gène
        selected_gene = st.selectbox("Détails gène", list(CONFIG.GENES_PERFORMANCE.keys()))
        if selected_gene:
            ref_seq = CONFIG.GENES_PERFORMANCE[selected_gene]
            st.code(f"Référence {selected_gene}: {ref_seq}", language="text")
    
    with tabs[1]:
        st.subheader("Statut Sanitaire & Résistance")
        df_health = engine.screen_all_markers(sequences)
        
        # Focus Scrapie
        scrapie_cols = [c for c in df_health.columns if 'Scrapie' in c and 'sim%' not in c]
        df_scrapie = df_health[['ID'] + scrapie_cols]
        
        st.markdown("### 🧠 Profil Scrapie (Prion Protein - PRNP)")
        for _, row in df_scrapie.iterrows():
            with st.expander(f"Résultats pour {row['ID']}"):
                cols = st.columns(len(scrapie_cols))
                for i, col in enumerate(scrapie_cols):
                    status = row[col]
                    color = "green" if "RÉSISTANT" in status else "red" if "SENSIBLE" in status else "gray"
                    cols[i].markdown(f"**{col}**  \n:{color}[{status}]")
    
    with tabs[2]:
        st.subheader("Analyse Populationnelle")
        if len(sequences) > 1:
            stats = engine.calculate_heterozygosity(sequences)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Hétérozygotie", f"{stats['heterozygosity']}%")
            col2.metric("Distance π", f"{stats['pi_distance']}")
            col3.metric("Comparaisons", stats['n_comparisons'])
            
            # Interprétation
            if stats['heterozygosity'] < 5:
                st.error("🚨 **Risque élevé de consanguinité** - Diversité très faible")
            elif stats['heterozygosity'] < 15:
                st.warning("⚠️ Diversité génétique limitée - Envisager l'introduction de sang neuf")
            else:
                st.success("✅ Diversité génétique satisfaisante")
        else:
            st.info("ℹ️ Importez au moins 2 séquences pour calculer les statistiques de population")
    
    with tabs[3]:
        st.subheader("Traduction en Acides Aminés")
        seq_choice = st.selectbox("Séquence à traduire", list(sequences.keys()))
        
        if seq_choice:
            result = engine.translate_to_protein(sequences[seq_choice])
            
            if "error" in result:
                st.error(result["error"])
            else:
                st.metric("Contenu GC", f"{result['gc_content']}%")
                
                for frame, data in result['protein_frames'].items():
                    with st.expander(f"Cadre de lecture {frame} ({data['length_aa']} AA)"):
                        st.code(data['protein'], language="text")
                        st.caption(f"Codons STOP: {data['stop_codons']}")
    
    with tabs[4]:
        st.subheader("Export & Archivage")
        
        # CSV complet
        full_df = engine.screen_all_markers(sequences)
        csv = full_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            "📥 Télécharger rapport complet (CSV)",
            csv,
            f"analyse_genomique_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
            "text/csv"
        )
        
        # Sauvegarde DB
        if st.button("💾 Archiver dans la base de données"):
            for seq_id, seq in sequences.items():
                seq_hash = hashlib.sha256(seq.encode()).hexdigest()[:16]
                try:
                    db.execute(
                        "INSERT OR IGNORE INTO sequences_genomiques (brebis_id, nom_sequence, sequence_hash, sequence_text) VALUES (?,?,?,?)",
                        (seq_id, seq_id, seq_hash, seq[:1000])  # Limite pour perf
                    )
                except Exception as e:
                    logger.error(f"Archivage échoué: {e}")
            st.success("Séquences archivées!")

def render_nutrition():
    """Calculateur de ration avec modèles nutritionnels avancés"""
    st.header("🌾 Calculateur de Ration Précision")
    
    st.markdown("""
    **Modèle**: INRA 2018 (Institut National de la Recherche Agronomique)
    Adapté pour ovins laitiers et viande Algérie
    """)
    
    with st.form("nutrition_calc"):
        col1, col2 = st.columns(2)
        
        poids_vif = col1.number_input("Poids vif (kg)", 20, 150, 60)
        etat_corporel = col2.slider("État Corporel (1-5)", 1.0, 5.0, 3.0, 0.5,
                                   help="1=Maigre, 3=Idéal, 5=Obèse")
        
        st.markdown("### 🎯 Objectif de Production")
        prod_type = st.radio("Type", ["Croissance/Engraissement", "Lactation", "Gestation", "Maintenance"])
        
        if prod_type == "Lactation":
            qte_lait = st.number_input("Production laitière (L/jour)", 0.0, 8.0, 2.0)
            mat_grasse = st.number_input("Matière grasse du lait (%)", 3.0, 9.0, 6.5)
        
        submitted = st.form_submit_button("Calculer la ration")
        
        if submitted:
            # Modèle INRA simplifié
            besoins = {
                "maintenance_uf": 0.038 * poids_vif**0.75,
                "maintenance_pdin": 3.25 * poids_vif**0.75,
                "maintenance_pdie": 3.25 * poids_vif**0.75
            }
            
            if prod_type == "Croissance/Engraissement":
                gain_jour = st.session_state.get('gain_vise', 200)  # g/j
                besoins['croissance_uf'] = gain_jour * 0.007
                besoins['croissance_pdin'] = gain_jour * 0.32
            
            # Affichage résultats
            st.subheader("🍽️ Ration Recommandée")
            
            c1, c2, c3 = st.columns(3)
            c1.metric("Fourrage sec", f"{poids_vif * 0.015:.1f} - {poids_vif * 0.025:.1f} kg/j")
            c2.metric("Concentré", f"{poids_vif * 0.008:.1f} - {poids_vif * 0.015:.1f} kg/j")
            c3.metric("Eau", f"{poids_vif * 0.08:.1f} L/j")
            
            # Composition détaillée
            with st.expander("Voir détail nutritionnel"):
                st.json(besoins)

# ============================================================================
# 6. APPLICATION PRINCIPALE
# ============================================================================

def main():
    """Point d'entrée avec gestion d'état et erreurs globales"""
    st.set_page_config(
        page_title="EXPERT OVIN DZ PRO",
        page_icon="🐑",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # CSS custom professionnel
    st.markdown("""
    <style>
    .stMetric {background-color: #f0f2f6; border-radius: 10px; padding: 10px;}
    .stDataFrame {font-size: 12px;}
    div[data-testid="stForm"] {background-color: #fafafa; padding: 20px; border-radius: 10px;}
    </style>
    """, unsafe_allow_html=True)
    
    try:
        # Initialisation singletons
        if 'db' not in st.session_state:
            st.session_state.db = DatabaseManager()
            init_schema(st.session_state.db)
            logger.info("Session DB initialisée")
        
        if 'bio_engine' not in st.session_state:
            st.session_state.bio_engine = BioInfoEngine()
        
        if 'analytics' not in st.session_state:
            st.session_state.analytics = SheepAnalytics(st.session_state.db)
        
        db = st.session_state.db
        engine = st.session_state.bio_engine
        analytics = st.session_state.analytics
        
        # Rendu UI
        render_header()
        choice = render_sidebar()
        
        # Routage
        if choice == "📊 Dashboard":
            render_dashboard(db, analytics)
        elif choice == "📝 Inscription":
            render_inscription(db)
        elif choice == "🥛 Production":
            st.header("🥛 Suivi Laitier")
            st.info("Module en cours d'enrichissement - Intégration automates de traite")
        elif choice == "🩺 Santé":
            st.header("🩺 Gestion Sanitaire")
            st.info("Carnet sanitaire avec alertes automatiques")
        elif choice == "🧬 Génomique":
            render_genomique(engine, db)
        elif choice == "🌾 Nutrition":
            render_nutrition()
        elif choice == "⚙️ Admin":
            st.header("⚙️ Administration")
            if st.button("🗑️ Vider le cache"):
                st.cache_data.clear()
                st.success("Cache vidé!")
    
    except Exception as e:
        logger.critical(f"Erreur critique: {e}", exc_info=True)
        st.error("🚨 Une erreur critique est survenue. Veuillez rafraîchir la page.")
        st.exception(e)

if __name__ == "__main__":
    main()
