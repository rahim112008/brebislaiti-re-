"""
EXPERT OVIN DZ PRO - VERSION INTEGRALE 2026.02
Système Tout-en-Un : Phénotypage, Lait, Génomique, Santé & Nutrition
Auteur : Gemini AI Collaborator
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sqlite3
import os
from datetime import datetime, date
from Bio import pairwise2
from Bio.Seq import Seq

# ============================================================================
# 1. MOTEUR DE BASE DE DONNÉES (PERSISTENCE)
# ============================================================================

class DatabaseManager:
    def __init__(self, db_path: str = "data/ovin_master_pro.db"):
        self.db_path = db_path
        if not os.path.exists('data'): os.makedirs('data')
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row

    def execute_query(self, query: str, params: tuple = ()):
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            st.error(f"Erreur SQL: {e}")
            return None

    def fetch_all_as_df(self, query: str, params: tuple = ()):
        return pd.read_sql_query(query, self.conn, params=params)

def init_database(db: DatabaseManager):
    """Initialise toutes les tables nécessaires au fonctionnement de l'app"""
    tables = [
        """CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT, race TEXT, poids REAL, note_mamelle INTEGER, 
            tour_poitrine REAL, longueur REAL, created_at DATE
        )""",
        """CREATE TABLE IF NOT EXISTS controle_laitier (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_controle DATE, quantite_lait REAL, 
            FOREIGN KEY (brebis_id) REFERENCES brebis (identifiant_unique)
        )""",
        """CREATE TABLE IF NOT EXISTS sante (
            id INTEGER PRIMARY KEY AUTOINCREMENT, brebis_id TEXT, 
            date_soin DATE, type_acte TEXT, produit TEXT, rappel_prevu DATE
        )"""
    ]
    for table_sql in tables: db.execute_query(table_sql)

# ============================================================================
# 2. MOTEUR BIOINFORMATIQUE & IA
# ============================================================================

class BioInfoEngine:
    # Références Génomiques (Performance & Santé)
    GENES_INTERET = {
        "FecB (Prolificité)": "GATGGTTCAAGTCCACAGTTTTA", 
        "MSTN (Muscle)": "AAGCTTGATTAGCAGGTTCCCGG",
        "CAST (Tendreté)": "TGGGGCCCAAGTCGATTGCAGAA",
        "DGAT1 (Lait)": "GCTAGCTAGCTAGCTGATCGATG"
    }
    GENES_PATHOLOGIES = {
        "Scrapie (Tremblante)": "TGGTACCCATAATCAGTGGAACA",
        "Arachnomélie": "CCGTAGCTAGCTGATCGATCGTA",
        "Hypotrichose": "TTAGCGCTAGCTAGCTAGCTAGC"
    }

    @staticmethod
    def filtrer_sequence(seq):
        """Nettoie la séquence des headers FASTA et espaces"""
        if ">" in seq: seq = "".join(seq.split('\n')[1:])
        return seq.upper().strip().replace(" ", "").replace("\r", "").replace("\n", "")

    @staticmethod
    def traduire_en_proteine(dna_seq):
        """Traduit la séquence ADN en Acides Aminés (Protéine)"""
        try:
            clean_dna = dna_seq[:(len(dna_seq)//3)*3]
            if not clean_dna: return "Séquence trop courte"
            return str(Seq(clean_dna).translate(to_stop=True))
        except: return "Erreur de traduction"

    @staticmethod
    def alignement_expert(seq_test, ref_seq):
        """Calcul de similarité par alignement local"""
        alignments = pairwise2.align.localxx(seq_test, ref_seq)
        if alignments:
            return round((alignments[0].score / len(ref_seq)) * 100, 2)
        return 0.0

# ============================================================================
# 3. INTERFACE UTILISATEUR PRINCIPALE
# ============================================================================

def main():
    st.set_page_config(page_title="EXPERT OVIN DZ PRO", layout="wide", page_icon="🐑")
    
    # Initialisation de la session et DB
    if 'db' not in st.session_state:
        st.session_state.db = DatabaseManager()
        init_database(st.session_state.db)
    
    db = st.session_state.db
    genomique = BioInfoEngine()

    # Sidebar Navigation
    st.sidebar.title("🐑 EXPERT OVIN DZ")
    st.sidebar.markdown("---")
    menu = [
        "📊 Dashboard Élite", 
        "📝 Inscription & Phénotype", 
        "📷 Scanner IA 1m", 
        "🥛 Contrôle Laitier", 
        "🩺 Santé & Vaccins", 
        "🧬 Génomique & NCBI", 
        "🌾 Nutrition Solo"
    ]
    choice = st.sidebar.radio("Navigation", menu)

    # --- 1. DASHBOARD ---
    if choice == "📊 Dashboard Élite":
        st.title("📊 Tableau de Bord des Performances")
        df_b = db.fetch_all_as_df("SELECT * FROM brebis")
        df_l = db.fetch_all_as_df("SELECT * FROM controle_laitier")
        
        if not df_b.empty:
            c1, c2, c3 = st.columns(3)
            c1.metric("Effectif Total", len(df_b))
            c2.metric("Poids Moyen", f"{round(df_b['poids'].mean(), 1)} kg")
            avg_lait = df_l['quantite_lait'].mean() if not df_l.empty else 0
            c3.metric("Moyenne Lait", f"{round(avg_lait, 2)} L")
            
            st.subheader("📋 Liste du Cheptel")
            st.dataframe(df_b, use_container_width=True)
            
            if not df_l.empty:
                st.subheader("📈 Évolution de la Production")
                fig = px.line(df_l, x='date_controle', y='quantite_lait', color='brebis_id', title="Courbe de Lactation")
                st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Aucun animal enregistré. Commencez par le module Inscription.")

    # --- 2. INSCRIPTION ---
    elif choice == "📝 Inscription & Phénotype":
        st.title("📝 Enregistrement Phénotypique")
        with st.form("form_inscription"):
            col1, col2 = st.columns(2)
            uid = col1.text_input("ID Boucle (Identifiant Unique)")
            nom = col1.text_input("Nom / Alias")
            race = col2.selectbox("Race", ["Ouled Djellal", "Rembi", "Hamra", "Lacaune", "Autre"])
            poids = col2.number_input("Poids (kg)", 10.0, 150.0, 50.0)
            
            st.markdown("🔍 **Mesures Morphométriques**")
            tp = st.number_input("Tour Poitrine (cm)", 40.0, 160.0, 85.0)
            lg = st.number_input("Longueur Corps (cm)", 30.0, 140.0, 75.0)
            note_m = st.slider("Note de Mamelle (1-10)", 1, 10, 5)
            
            if st.form_submit_button("Sauvegarder l'animal"):
                db.execute_query(
                    "INSERT INTO brebis (identifiant_unique, nom, race, poids, note_mamelle, tour_poitrine, longueur, created_at) VALUES (?,?,?,?,?,?,?,?)",
                    (uid, nom, race, poids, note_m, tp, lg, date.today())
                )
                st.success(f"L'animal {uid} a été ajouté avec succès.")

    # --- 3. SCANNER IA ---
    elif choice == "📷 Scanner IA 1m":
        st.title("📸 Scanner IA Morphométrique")
        st.info("Prenez une photo latérale de l'animal avec un étalon de 1 mètre placé à côté.")
        st.camera_input("Scanner")
        st.warning("Module de calcul de pixels en attente de calibration avec l'étalon.")

    # --- 4. LAIT ---
    elif choice == "🥛 Contrôle Laitier":
        st.title("🥛 Suivi de Production Laitière")
        with st.form("form_lait"):
            id_lait = st.text_input("Scanner l'ID de la brebis")
            qte_l = st.number_input("Quantité de lait (L)", 0.0, 12.0, 1.5)
            dt_l = st.date_input("Date du contrôle", date.today())
            if st.form_submit_button("Valider la traite"):
                db.execute_query("INSERT INTO controle_laitier (brebis_id, date_controle, quantite_lait) VALUES (?,?,?)",
                                 (id_lait, dt_l, qte_l))
                st.success("Donnée enregistrée.")

    # --- 5. SANTÉ ---
    elif choice == "🩺 Santé & Vaccins":
        st.title("🩺 Gestion Sanitaire")
        with st.expander("➕ Enregistrer un nouvel acte"):
            with st.form("form_sante"):
                id_s = st.text_input("ID de l'animal")
                type_a = st.selectbox("Acte", ["Vaccination", "Déparasitage", "Traitement Curatif"])
                prod = st.text_input("Médicament / Produit utilisé")
                rappel = st.date_input("Date de rappel (si applicable)", date.today() + timedelta(days=30))
                if st.form_submit_button("Ajouter au carnet"):
                    db.execute_query("INSERT INTO sante (brebis_id, date_soin, type_acte, produit, rappel_prevu) VALUES (?,?,?,?,?)",
                                     (id_s, date.today(), type_a, prod, rappel))
                    st.success("Soin enregistré.")
        
        st.subheader("📅 Historique des Soins")
        df_sante = db.fetch_all_as_df("SELECT * FROM sante")
        st.table(df_sante)

   # --- 6. GÉNOMIQUE ---
    elif choice == "🧬 Génomique & NCBI":
        st.title("🧬 Laboratoire de Génomique & Bio-informatique")
        
        # Zone de saisie flexible (Simple ou Multi-FASTA)
        dna_txt = st.text_area("Collez vos séquences ADN (Format FASTA, Multi-FASTA ou Brut)", height=200)
        
        if dna_txt:
            # 1. Préparation des données
            if dna_txt.count(">") > 1:
                # Mode Multi-Fasta pour criblage de masse
                data_dict = genomique.extraire_multi_fasta(dna_txt)
                is_multi = True
            else:
                # Mode individuel
                seq_val = genomique.filtrer_sequence(dna_txt)
                data_dict = {"Individu_Unique": seq_val}
                is_multi = False

            # 2. Création des onglets enrichis
            t_perf, t_patho, t_pop, t_trad = st.tabs([
                "🎯 Performance SNP", 
                "⚠️ Santé & Résistance ARR", 
                "📊 Diversité Élevage", 
                "🔬 Analyse Moléculaire"
            ])
            
            # --- ONGLET PERFORMANCE ---
            with t_perf:
                st.subheader("Criblage des Marqueurs de Production")
                results_perf = []
                for name, sequence in data_dict.items():
                    row = {"ID": name}
                    for gene, ref in genomique.GENES_INTERET.items():
                        score = genomique.alignement_expert(sequence, ref)
                        status = "OUI" if score > 85 else "NON"
                        row[gene] = f"{status} ({score}%)"
                    results_perf.append(row)
                
                df_perf = pd.DataFrame(results_perf)
                st.dataframe(df_perf, use_container_width=True)
                
                # Export pour l'éleveur
                csv_perf = df_perf.to_csv(index=False).encode('utf-8')
                st.download_button("📥 Télécharger le Catalogue Performance (CSV)", csv_perf, "performance_elevage.csv", "text/csv")

            # --- ONGLET SANTÉ & RÉSISTANCE ---
            with t_patho:
                st.subheader("🛡️ Statut Sanitaire & Résistance Tremblante")
                results_sante = []
                for name, sequence in data_dict.items():
                    row = {"ID": name}
                    for path, ref in genomique.GENES_SANTE.items():
                        score = genomique.alignement_expert(sequence, ref)
                        if score > 85:
                            res = "POSITIF" if "VRQ" in path or "Scrapie" not in path else "RÉSISTANT (ARR)"
                        else:
                            res = "NÉGATIF"
                        row[path] = res
                    results_sante.append(row)
                
                df_sante = pd.DataFrame(results_sante)
                st.table(df_sante)
                
                

            # --- ONGLET DIVERSITÉ (LE NOUVEAU) ---
            with t_pop:
                st.subheader("📊 Étude de Population & Consanguinité")
                if not is_multi:
                    st.info("ℹ️ Pour calculer le taux d'hétérozygotie, veuillez coller au moins 2 séquences au format Multi-FASTA.")
                else:
                    score_h = genomique.calculer_heterozygotie(data_dict)
                    col1, col2 = st.columns(2)
                    col1.metric("Indice de Diversité (Hétérozygotie)", f"{score_h}%")
                    
                    if score_h < 10:
                        col2.error("⚠️ Risque de consanguinité élevé dans cet élevage.")
                    else:
                        col2.success("✅ Bonne variabilité génétique détectée.")
                    
                    # Graphique de proximité
                    
                    st.write("**Note :** Plus l'hétérozygotie est élevée, plus l'élevage est résistant aux maladies et productif sur le long terme.")

            # --- ONGLET TRADUCTION ---
            with t_trad:
                st.subheader("Séquençage Protéique")
                # On prend le premier de la liste pour la traduction
                premier_id = list(data_dict.keys())[0]
                st.write(f"Traduction de l'individu : **{premier_id}**")
                st.code(genomique.traduire_en_proteine(data_dict[premier_id]), language="text")
                

    # --- 7. NUTRITION ---
    elif choice == "🌾 Nutrition Solo":
        st.title("🌾 Calculateur de Ration de Précision")
        p_indiv = st.number_input("Poids de l'animal (kg)", 10, 150, 60)
        c1, c2 = st.columns(2)
        c1.write(f"🌾 **Concentré (Orge) :** {p_indiv * 0.012:.2f} kg/jour")
        c2.write(f"🌿 **Fourrage (Foin/Luzerne) :** {p_indiv * 0.02:.2f} kg/jour")

if __name__ == "__main__":
    from datetime import timedelta
    main()
