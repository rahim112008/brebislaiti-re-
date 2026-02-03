"""
OVIN MANAGER PRO - Version Streamlit Cloud Compatible
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import sqlite3
import numpy as np
import json
import io

# ========== CONFIGURATION ==========

st.set_page_config(
    page_title="Ovin Manager Pro",
    page_icon="🐑",
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
    .dna-sequence {
        font-family: 'Courier New', monospace;
        background-color: #f0f0f0;
        padding: 10px;
        border-radius: 5px;
        margin: 5px 0;
        font-size: 0.9em;
        letter-spacing: 1px;
    }
</style>
""", unsafe_allow_html=True)

# ========== INITIALISATION BASE DE DONNÉES ==========

def init_database():
    """Initialise la base de données SQLite"""
    conn = sqlite3.connect('ovin_manager.db', check_same_thread=False)
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
    
    # Table données génomiques simplifiée
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS donnees_genomiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            gene_nom TEXT,
            sequence_adn TEXT,
            chromosome TEXT,
            genotype TEXT,
            date_analyse DATE,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
    ''')
    
    conn.commit()
    return conn

# Connexion à la base de données
conn = init_database()

# ========== MODULES DE FONCTIONS ==========

class GeneticAnalyzer:
    """Analyseur génétique simplifié"""
    
    @staticmethod
    def analyze_sequence(sequence):
        """Analyse une séquence ADN"""
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        
        if len(seq) == 0:
            return {"error": "Séquence vide"}
        
        # Composition
        composition = {
            'A': seq.count('A'),
            'T': seq.count('T'),
            'C': seq.count('C'),
            'G': seq.count('G'),
            'N': seq.count('N') + seq.count('X')
        }
        
        total = len(seq)
        gc_content = (composition['G'] + composition['C']) / total * 100 if total > 0 else 0
        
        # Motifs
        motifs = {
            'start_codon': seq.count('ATG'),
            'stop_codons': seq.count('TAA') + seq.count('TAG') + seq.count('TGA'),
            'cpgi_sites': seq.count('CG')
        }
        
        return {
            'length': total,
            'composition': composition,
            'gc_content': round(gc_content, 2),
            'at_content': round(100 - gc_content, 2),
            'motifs': motifs,
            'gc_skew': GeneticAnalyzer.calculate_gc_skew(seq)
        }
    
    @staticmethod
    def calculate_gc_skew(seq):
        """Calcule le GC skew"""
        g_count = seq.count('G')
        c_count = seq.count('C')
        total = g_count + c_count
        return round((g_count - c_count) / total, 3) if total > 0 else 0
    
    @staticmethod
    def translate_dna_to_protein(dna_seq):
        """Traduction ADN -> Protéine (simplifiée)"""
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
        for i in range(0, len(dna_seq)-2, 3):
            codon = dna_seq[i:i+3]
            protein += codon_table.get(codon, 'X')
        
        return protein

# ========== PAGES DE L'APPLICATION ==========

def page_accueil():
    """Page d'accueil"""
    st.markdown('<h1 class="main-header">🐑 Ovin Manager Pro</h1>', unsafe_allow_html=True)
    st.markdown("*Application scientifique de gestion et d'analyse d'élevage ovin laitier*")
    
    # Métriques rapides
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM brebis")
        total = cursor.fetchone()[0]
        st.metric("Total Brebis", total)
    
    with col2:
        cursor.execute("SELECT COUNT(*) FROM brebis WHERE sexe = 'F'")
        femelles = cursor.fetchone()[0]
        st.metric("Brebis Femelles", femelles)
    
    with col3:
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        gestations = cursor.fetchone()[0]
        st.metric("Gestations", gestations)
    
    with col4:
        cursor.execute("SELECT AVG(quantite_litre) FROM production_lait")
        prod = cursor.fetchone()[0] or 0
        st.metric("Lait moyen/j", f"{prod:.1f}L")
    
    # Graphiques
    st.markdown("### 📊 Vue d'ensemble")
    
    col_graph1, col_graph2 = st.columns(2)
    
    with col_graph1:
        cursor.execute("SELECT race, COUNT(*) as count FROM brebis GROUP BY race")
        data_race = cursor.fetchall()
        
        if data_race:
            df_race = pd.DataFrame(data_race, columns=['Race', 'Nombre'])
            fig = px.pie(df_race, values='Nombre', names='Race', 
                        title="Répartition par race", hole=0.4)
            st.plotly_chart(fig, use_container_width=True)
    
    with col_graph2:
        cursor.execute("SELECT statut, COUNT(*) as count FROM brebis GROUP BY statut")
        data_statut = cursor.fetchall()
        
        if data_statut:
            df_statut = pd.DataFrame(data_statut, columns=['Statut', 'Nombre'])
            fig = px.bar(df_statut, x='Statut', y='Nombre', 
                        title="Répartition par statut", color='Statut')
            st.plotly_chart(fig, use_container_width=True)

def page_gestion_brebis():
    """Gestion des brebis"""
    st.markdown('<h2 class="section-header">📊 Gestion des Brebis</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["📋 Liste", "➕ Ajouter", "🔍 Rechercher", "📤 Exporter"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM brebis ORDER BY id DESC")
        brebis_data = cursor.fetchall()
        
        if brebis_data:
            columns = [desc[0] for desc in cursor.description]
            df = pd.DataFrame(brebis_data, columns=columns)
            
            # Recherche
            recherche = st.text_input("🔍 Rechercher une brebis:")
            if recherche:
                mask = df.apply(lambda row: row.astype(str).str.contains(recherche, case=False).any(), axis=1)
                df = df[mask]
            
            st.dataframe(df[['identifiant_unique', 'nom', 'race', 'sexe', 'statut', 'poids', 'date_naissance']], 
                        use_container_width=True)
        else:
            st.info("Aucune brebis enregistrée")
    
    with tab2:
        with st.form("form_ajout_brebis", clear_on_submit=True):
            col1, col2 = st.columns(2)
            
            with col1:
                identifiant = st.text_input("Identifiant Unique*", placeholder="BRB-2024-001")
                nom = st.text_input("Nom", placeholder="Bella")
                date_naissance = st.date_input("Date naissance", value=date.today() - timedelta(days=365))
                race = st.selectbox("Race", ["Ouled Djellal", "Razè", "Hamra", "D'man", "Saharienne", "Croisé", "Autre"])
            
            with col2:
                sexe = st.radio("Sexe", ["F", "M"], horizontal=True)
                statut = st.selectbox("Statut", ["active", "retrait", "malade", "vendu"])
                poids = st.number_input("Poids (kg)", min_value=0.0, value=50.0)
                notes = st.text_area("Notes")
            
            if st.form_submit_button("✅ Ajouter la brebis", type="primary"):
                if identifiant:
                    try:
                        cursor = conn.cursor()
                        cursor.execute('''
                            INSERT INTO brebis (identifiant_unique, nom, date_naissance, race, sexe, statut, poids, notes)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (identifiant, nom, date_naissance.isoformat(), race, sexe, statut, poids, notes))
                        conn.commit()
                        st.success(f"✅ Brebis {nom} ajoutée avec succès!")
                        st.balloons()
                    except Exception as e:
                        st.error(f"Erreur: {e}")
                else:
                    st.warning("L'identifiant unique est obligatoire")
    
    with tab3:
        st.info("Fonction de recherche avancée en développement...")
    
    with tab4:
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM brebis")
        data = cursor.fetchall()
        
        if data:
            columns = [desc[0] for desc in cursor.description]
            df = pd.DataFrame(data, columns=columns)
            
            st.write("Exporter les données au format:")
            
            col_exp1, col_exp2 = st.columns(2)
            
            with col_exp1:
                if st.button("📄 CSV"):
                    csv = df.to_csv(index=False)
                    st.download_button(
                        label="Télécharger CSV",
                        data=csv,
                        file_name="brebis_export.csv",
                        mime="text/csv"
                    )
            
            with col_exp2:
                if st.button("📋 JSON"):
                    json_data = df.to_json(orient='records', indent=2)
                    st.download_button(
                        label="Télécharger JSON",
                        data=json_data,
                        file_name="brebis_export.json",
                        mime="application/json"
                    )

def page_genetique():
    """Module génétique"""
    st.markdown('<h2 class="section-header">🧬 Module Génétique</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4 = st.tabs(["🧬 Analyse ADN", "🔍 Gènes Ovin", "📊 SNP", "💾 Données"])
    
    with tab1:
        st.markdown("### Analyse de séquences ADN")
        
        col_seq1, col_seq2 = st.columns([2, 1])
        
        with col_seq1:
            sequence = st.text_area(
                "Collez votre séquence ADN:",
                height=200,
                placeholder="ATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGC..."
            )
            
            if st.button("🔬 Analyser la séquence", type="primary"):
                if sequence:
                    result = GeneticAnalyzer.analyze_sequence(sequence)
                    
                    if "error" not in result:
                        st.success(f"✅ Séquence analysée: {result['length']} bp")
                        
                        # Métriques
                        col_met1, col_met2, col_met3, col_met4 = st.columns(4)
                        
                        with col_met1:
                            st.metric("Longueur", f"{result['length']} bp")
                        with col_met2:
                            st.metric("% GC", f"{result['gc_content']}%")
                        with col_met3:
                            st.metric("Codons START", result['motifs']['start_codon'])
                        with col_met4:
                            st.metric("GC Skew", f"{result['gc_skew']:.3f}")
                        
                        # Graphique composition
                        fig = px.pie(
                            values=list(result['composition'].values()),
                            names=list(result['composition'].keys()),
                            title="Composition nucléotidique"
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Traduction protéique
                        if result['length'] >= 3:
                            protein = GeneticAnalyzer.translate_dna_to_protein(sequence[:300])  # Premiers 300 bp
                            st.markdown("**Traduction protéique (premiers 100 acides aminés):**")
                            st.code(protein[:100])
                    else:
                        st.error(result["error"])
        
        with col_seq2:
            st.markdown("#### Exemples")
            examples = {
                "Court (100bp)": "ATCG" * 25,
                "Séquence MSTN": "ATGGCCATTGAACAGAAACCAACCTACCCCGAGAACAGCTTTGAGGACAGCCTGGGCCGCATGG",
                "Riche en GC": "GGCCGGCC" * 20
            }
            
            for name, seq in examples.items():
                if st.button(f"Charger: {name}"):
                    st.session_state.example_seq = seq
                    st.rerun()
            
            if 'example_seq' in st.session_state:
                sequence = st.session_state.example_seq
    
    with tab2:
        st.markdown("### Base de données des gènes ovins")
        
        genes_db = {
            "MSTN": {
                "nom": "Myostatine",
                "chromosome": "2",
                "fonction": "Régulateur de croissance musculaire",
                "phénotype": "Double-muscling",
                "mutations": ["g.6723G>A", "c.939G>A"]
            },
            "PRNP": {
                "nom": "Protéine Prion",
                "chromosome": "13",
                "fonction": "Résistance à la tremblante",
                "phénotype": "Résistance aux ESST",
                "mutations": ["codon 136", "codon 154", "codon 171"]
            },
            "DGAT1": {
                "nom": "Diacylglycérol acyltransférase",
                "chromosome": "14",
                "fonction": "Synthèse des triglycérides",
                "phénotype": "Teneur en matière grasse du lait",
                "mutations": ["K232A"]
            },
            "GDF9": {
                "nom": "Growth Differentiation Factor 9",
                "chromosome": "5",
                "fonction": "Fertilité femelle",
                "phénotype": "Prolificité",
                "mutations": ["G1", "G4"]
            }
        }
        
        gene_selected = st.selectbox("Sélectionnez un gène", list(genes_db.keys()))
        
        if gene_selected:
            gene_info = genes_db[gene_selected]
            
            col_gene1, col_gene2 = st.columns(2)
            
            with col_gene1:
                st.markdown(f"#### {gene_selected} - {gene_info['nom']}")
                st.write(f"**Chromosome:** {gene_info['chromosome']}")
                st.write(f"**Fonction:** {gene_info['fonction']}")
                st.write(f"**Phénotype associé:** {gene_info['phénotype']}")
            
            with col_gene2:
                st.markdown("#### Mutations connues")
                for mut in gene_info['mutations']:
                    st.write(f"- {mut}")
                
                if st.button(f"📥 Importer {gene_selected}"):
                    # Stocker dans la base
                    cursor = conn.cursor()
                    cursor.execute('''
                        INSERT INTO donnees_genomiques (gene_nom, chromosome, date_analyse)
                        VALUES (?, ?, ?)
                    ''', (gene_selected, gene_info['chromosome'], date.today().isoformat()))
                    conn.commit()
                    st.success(f"Gène {gene_selected} importé!")
    
    with tab3:
        st.markdown("### Analyse des marqueurs SNP")
        
        # Données simulées
        snp_data = pd.DataFrame({
            'SNP': ['rs123456', 'rs789012', 'rs345678', 'rs901234'],
            'Chromosome': ['2', '6', '14', '19'],
            'Position': [123456, 789012, 345678, 901234],
            'Allèle majeur': ['A', 'G', 'C', 'T'],
            'Allèle mineur': ['G', 'A', 'T', 'C'],
            'MAF': [0.42, 0.18, 0.33, 0.25],
            'Gène': ['MSTN', 'PRNP', 'DGAT1', 'GDF9']
        })
        
        st.dataframe(snp_data)
        
        # Graphique des fréquences
        fig = px.bar(snp_data, x='SNP', y='MAF', color='Gène',
                    title="Fréquence de l'allèle mineur (MAF) par SNP")
        st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        st.markdown("### Gestion des données génomiques")
        
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM donnees_genomiques")
        count = cursor.fetchone()[0]
        
        st.metric("Séquences stockées", count)
        
        if count > 0:
            cursor.execute("SELECT * FROM donnees_genomiques ORDER BY id DESC LIMIT 10")
            data = cursor.fetchall()
            columns = [desc[0] for desc in cursor.description]
            df = pd.DataFrame(data, columns=columns)
            st.dataframe(df)

def page_analyse_lait():
    """Analyse laitière"""
    st.markdown('<h2 class="section-header">🥛 Analyse Laitière</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📝 Saisie", "📈 Graphiques", "📊 Statistiques"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom FROM brebis WHERE sexe = 'F'")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            with st.form("form_production_lait"):
                brebis_selected = st.selectbox(
                    "Sélectionner une brebis",
                    [f"{b[1]} (ID: {b[0]})" for b in brebis_list]
                )
                
                date_mesure = st.date_input("Date de mesure", value=date.today())
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    quantite = st.number_input("Quantité (L)", 0.0, 10.0, 2.5, 0.1)
                
                with col2:
                    mg = st.number_input("Matière grasse %", 0.0, 20.0, 6.5, 0.1)
                
                with col3:
                    proteine = st.number_input("Protéine %", 0.0, 20.0, 5.2, 0.1)
                
                notes = st.text_area("Notes")
                
                if st.form_submit_button("💾 Enregistrer", type="primary"):
                    brebis_id = int(brebis_selected.split("ID: ")[1].rstrip(")"))
                    
                    cursor.execute('''
                        INSERT INTO production_lait (brebis_id, date_mesure, quantite_litre, 
                                                   taux_matiere_grasse, taux_proteine, notes)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (brebis_id, date_mesure.isoformat(), quantite, mg, proteine, notes))
                    conn.commit()
                    st.success("✅ Production enregistrée!")
        else:
            st.warning("Aucune brebis femelle enregistrée")
    
    with tab2:
        cursor.execute("""
            SELECT date_mesure, AVG(quantite_litre) as lait, 
                   AVG(taux_matiere_grasse) as mg, AVG(taux_proteine) as proteine
            FROM production_lait
            WHERE date_mesure > date('now', '-30 days')
            GROUP BY date_mesure
            ORDER BY date_mesure
        """)
        
        data = cursor.fetchall()
        
        if data:
            df = pd.DataFrame(data, columns=['Date', 'Lait (L)', 'MG (%)', 'Protéine (%)'])
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=df['Date'], y=df['Lait (L)'], name='Lait (L)'))
            fig.add_trace(go.Scatter(x=df['Date'], y=df['MG (%)'], name='MG (%)', yaxis='y2'))
            
            fig.update_layout(
                title="Production laitière - 30 derniers jours",
                yaxis=dict(title="Lait (L)"),
                yaxis2=dict(title="%", overlaying='y', side='right'),
                hovermode='x unified'
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        cursor.execute("""
            SELECT b.race, 
                   AVG(p.quantite_litre) as lait_moyen,
                   AVG(p.taux_matiere_grasse) as mg_moyen,
                   AVG(p.taux_proteine) as proteine_moyenne,
                   COUNT(*) as nb_mesures
            FROM production_lait p
            JOIN brebis b ON p.brebis_id = b.id
            GROUP BY b.race
            HAVING nb_mesures >= 3
        """)
        
        data_stats = cursor.fetchall()
        
        if data_stats:
            df_stats = pd.DataFrame(data_stats, 
                                   columns=['Race', 'Lait moyen (L)', 'MG moyenne (%)', 
                                           'Protéine moyenne (%)', 'Mesures'])
            
            st.dataframe(df_stats)
            
            # Graphique radar
            fig = go.Figure()
            
            for idx, row in df_stats.iterrows():
                fig.add_trace(go.Scatterpolar(
                    r=[row['Lait moyen (L)'], row['MG moyenne (%)'], row['Protéine moyenne (%)']],
                    theta=['Lait', 'MG', 'Protéine'],
                    fill='toself',
                    name=row['Race']
                ))
            
            fig.update_layout(
                polar=dict(radialaxis=dict(visible=True)),
                title="Comparaison par race"
            )
            
            st.plotly_chart(fig, use_container_width=True)

def page_gestation():
    """Suivi des gestations"""
    st.markdown('<h2 class="section-header">🤰 Suivi des Gestations</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["📅 Calendrier", "➕ Nouvelle", "📊 Statistiques"])
    
    with tab1:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT g.*, b.nom, b.race, b.identifiant_unique
            FROM gestations g 
            JOIN brebis b ON g.brebis_id = b.id 
            WHERE g.statut = 'en_cours'
            ORDER BY g.date_mise_bas_prevu
        ''')
        
        gestations = cursor.fetchall()
        
        if gestations:
            cols = ['id', 'brebis_id', 'date_eponge', 'date_mise_bas_prevu', 
                   'nombre_agneaux_prevus', 'statut', 'nom', 'race', 'identifiant_unique']
            
            today = date.today()
            
            for gest in gestations:
                gest_dict = dict(zip(cols, gest))
                mise_bas = datetime.strptime(gest_dict['date_mise_bas_prevu'], '%Y-%m-%d').date()
                jours_restants = (mise_bas - today).days
                
                # Couleur selon proximité
                if jours_restants < 0:
                    color = "🔴"
                elif jours_restants <= 7:
                    color = "🟠"
                elif jours_restants <= 30:
                    color = "🟡"
                else:
                    color = "🟢"
                
                with st.expander(f"{color} {gest_dict['nom']} - {gest_dict['identifiant_unique']}"):
                    col_info1, col_info2 = st.columns(2)
                    
                    with col_info1:
                        st.write(f"**Race:** {gest_dict['race']}")
                        st.write(f"**Date éponge:** {gest_dict['date_eponge']}")
                        st.write(f"**Jours gestation:** {(today - datetime.strptime(gest_dict['date_eponge'], '%Y-%m-%d').date()).days}")
                    
                    with col_info2:
                        st.write(f"**Mise bas prévue:** {gest_dict['date_mise_bas_prevu']}")
                        st.write(f"**Jours restants:** {jours_restants}")
                        st.write(f"**Agneaux prévus:** {gest_dict['nombre_agneaux_prevus']}")
                    
                    # Progression
                    progression = min(1.0, (today - datetime.strptime(gest_dict['date_eponge'], '%Y-%m-%d').date()).days / 150)
                    st.progress(progression)
                    
                    if st.button(f"✅ Mise bas réalisée", key=f"done_{gest_dict['id']}"):
                        cursor.execute("UPDATE gestations SET statut = 'termine' WHERE id = ?", 
                                     (gest_dict['id'],))
                        conn.commit()
                        st.success("Statut mis à jour!")
                        st.rerun()
        else:
            st.info("Aucune gestation en cours")
    
    with tab2:
        cursor.execute("SELECT id, nom FROM brebis WHERE sexe = 'F' AND statut = 'active'")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            with st.form("form_nouvelle_gestation"):
                brebis_selected = st.selectbox("Brebis", 
                                              [f"{b[1]} (ID: {b[0]})" for b in brebis_list])
                date_eponge = st.date_input("Date d'éponge", value=date.today())
                nb_agneaux = st.number_input("Nombre d'agneaux prévus", 1, 4, 1)
                notes = st.text_area("Notes")
                
                if st.form_submit_button("📅 Enregistrer", type="primary"):
                    brebis_id = int(brebis_selected.split("ID: ")[1].rstrip(")"))
                    date_mise_bas = date_eponge + timedelta(days=150)
                    
                    cursor.execute('''
                        INSERT INTO gestations (brebis_id, date_eponge, date_mise_bas_prevu, 
                                              nombre_agneaux_prevus, notes)
                        VALUES (?, ?, ?, ?, ?)
                    ''', (brebis_id, date_eponge.isoformat(), 
                         date_mise_bas.isoformat(), nb_agneaux, notes))
                    conn.commit()
                    st.success("✅ Gestation enregistrée!")
        else:
            st.warning("Aucune brebis femelle active disponible")
    
    with tab3:
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
        en_cours = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'termine'")
        terminees = cursor.fetchone()[0]
        
        col_stat1, col_stat2, col_stat3 = st.columns(3)
        
        with col_stat1:
            st.metric("En cours", en_cours)
        with col_stat2:
            st.metric("Terminées", terminees)
        with col_stat3:
            total = en_cours + terminees
            taux = (terminees / total * 100) if total > 0 else 0
            st.metric("Taux réussite", f"{taux:.1f}%")

def page_parametres():
    """Paramètres"""
    st.markdown('<h2 class="section-header">⚙️ Paramètres</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["Général", "Base de données", "À propos"])
    
    with tab1:
        st.markdown("### Paramètres généraux")
        
        col1, col2 = st.columns(2)
        
        with col1:
            langue = st.selectbox("Langue", ["Français", "English", "Español"])
            unite_poids = st.radio("Unité de poids", ["kg", "lbs"])
            format_date = st.selectbox("Format date", ["JJ/MM/AAAA", "AAAA-MM-JJ"])
        
        with col2:
            notifications = st.checkbox("Activer notifications", True)
            theme = st.selectbox("Thème", ["Clair", "Sombre", "Auto"])
            auto_save = st.checkbox("Sauvegarde auto", True)
        
        if st.button("💾 Sauvegarder"):
            st.success("Paramètres sauvegardés!")
    
    with tab2:
        st.markdown("### Gestion base de données")
        
        col_db1, col_db2, col_db3 = st.columns(3)
        
        with col_db1:
            if st.button("🗑️ Vider cache"):
                st.warning("Cette action est irréversible!")
        
        with col_db2:
            if st.button("💾 Sauvegarde"):
                # Exporter données
                cursor = conn.cursor()
                tables = ['brebis', 'gestations', 'production_lait', 'donnees_genomiques']
                backup = {}
                
                for table in tables:
                    cursor.execute(f"SELECT * FROM {table}")
                    backup[table] = cursor.fetchall()
                
                st.download_button(
                    label="📥 Télécharger backup",
                    data=json.dumps(backup, indent=2),
                    file_name=f"backup_ovin_{date.today()}.json",
                    mime="application/json"
                )
        
        with col_db3:
            uploaded_file = st.file_uploader("Restaurer backup", type=['json'])
            if uploaded_file:
                st.warning("La restauration écrasera les données actuelles!")
    
    with tab3:
        st.markdown("### À propos")
        
        st.markdown("""
        **Ovin Manager Pro** v2.0
        
        Application de gestion scientifique d'élevage ovin laitier
        
        **Fonctionnalités:**
        - Gestion complète du troupeau
        - Analyse génétique
        - Suivi laitière
        - Gestion des gestations
        - Statistiques avancées
        
        **Technologies:**
        - Python 3.10+
        - Streamlit
        - SQLite
        - Plotly
        
        **Développeur:** Rahim112008
        **Licence:** MIT
        """)

# ========== NAVIGATION PRINCIPALE ==========

# Titre principal (toujours affiché)
st.markdown('<h1 class="main-header">🐑 Ovin Manager Pro</h1>', unsafe_allow_html=True)

# Sidebar navigation
with st.sidebar:
    st.markdown("### 📍 Navigation")
    
    page = st.radio(
        "Menu Principal",
        ["🏠 Accueil", 
         "📊 Gestion Brebis", 
         "🧬 Génétique",
         "🥛 Analyse Lait",
         "🤰 Gestations",
         "⚙️ Paramètres"]
    )
    
    st.markdown("---")
    st.markdown("### 📊 Statistiques")
    
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM brebis")
    total = cursor.fetchone()[0]
    st.metric("Brebis", total)
    
    cursor.execute("SELECT COUNT(*) FROM gestations WHERE statut = 'en_cours'")
    gest = cursor.fetchone()[0]
    st.metric("Gestations", gest)

# Affichage de la page sélectionnée
if page == "🏠 Accueil":
    page_accueil()
elif page == "📊 Gestion Brebis":
    page_gestion_brebis()
elif page == "🧬 Génétique":
    page_genetique()
elif page == "🥛 Analyse Lait":
    page_analyse_lait()
elif page == "🤰 Gestations":
    page_gestation()
elif page == "⚙️ Paramètres":
    page_parametres()

# Pied de page
st.markdown("---")
st.caption(f"🐑 Ovin Manager Pro v2.0 | {date.today()} | © 2024")
