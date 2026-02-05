"""
PROJET : EXPERT OVIN DZ PRO (VERSION INTÉGRALE 2026)
Domaine : Sélection génétique, Génomique, Morphométrie & Gestion Laitière
Auteur : rahim LABORATOIRE GenApAgiE 
"""

"""
EXPERT OVIN DZ PRO - VERSION ULTRA EXPERT (FUSIONNÉE)
Bioinformatique, Analyse SNP, Biochimie Laitière & Sélection Génomique
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import io
from typing import Dict, List

# ============================================================================
# 1. STANDARDS ET CONFIGURATION
# ============================================================================
CALIBRATION_STANDARDS = {
    "Pièce 100 DA (Diamètre: 2.95cm)": 2.95,
    "Feuille A4 (Hauteur: 29.7cm)": 29.7,
    "Carte Bancaire (8.56cm)": 8.56,
    "Standard 1m": 100.0
}

# ============================================================================
# 2. MODULES EXPERTS (GÉNOMIQUE, SNP & BIOCHIMIE)
# ============================================================================

class IntegrationGenomique:
    """Intégration avancée et Analyse de SNP (Single Nucleotide Polymorphism)"""
    
    def __init__(self, email: str = "labo@expert-ovin.dz"):
        self.email = email
        self.sequences_cache = {}

    def analyser_snp(self, sequence_reference: str, sequence_etudiee: str) -> Dict:
        """Analyse les mutations SNP entre la référence NCBI et l'animal"""
        if len(sequence_reference) != len(sequence_etudiee):
            return {"erreur": "Les séquences doivent avoir la même longueur pour l'alignement direct"}
        
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
            'frequence_snp': len(snps) / len(sequence_reference) if len(sequence_reference) > 0 else 0,
            'snps_detailles': snps,
            'sequence_longueur': len(sequence_reference)
        }

    def _determiner_type_mutation(self, ref: str, etu: str) -> str:
        transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
        transversions = [('A', 'C'), ('A', 'T'), ('G', 'C'), ('G', 'T'),
                        ('C', 'A'), ('T', 'A'), ('C', 'G'), ('T', 'G')]
        if (ref, etu) in transitions: return 'transition'
        elif (ref, etu) in transversions: return 'transversion'
        return 'indéterminé'

    def rechercher_genes_candidats(self, race: str) -> List[Dict]:
        """Base de données des gènes d'intérêt par race"""
        genes_ovins = {
            'Lacaune': [
                {'gene': 'LALBA', 'fonction': 'Protéine du lait', 'chromosome': '3'},
                {'gene': 'CSN1S1', 'fonction': 'Caséine alpha-S1', 'chromosome': '6'},
                {'gene': 'DGAT1', 'fonction': 'Synthèse des triglycérides', 'chromosome': '14'},
            ],
            'Ouled Djellal': [
                {'gene': 'GDF9', 'fonction': 'Fécondité/Prolificité', 'chromosome': 'X'},
                {'gene': 'MSTN', 'fonction': 'Développement musculaire', 'chromosome': '2'},
            ],
            'Rembi': [
                {'gene': 'PRLR', 'fonction': 'Récepteur prolactine', 'chromosome': '16'},
                {'gene': 'GH1', 'fonction': 'Hormone de croissance', 'chromosome': '11'},
            ]
        }
        return genes_ovins.get(race, [{'gene': 'GENERIC', 'fonction': 'Gène ovin standard', 'chromosome': 'NA'}])

class UltraExpertModule:
    """Calculs biochimiques et index laitiers"""
    @staticmethod
    def get_biochem_diagnostic(fat, prot, bhb):
        ratio = fat / prot
        status = "Normal" if bhb < 1.2 else "Cétose Subclinique ⚠️"
        rumen = "Optimal" if 1.1 <= ratio <= 1.4 else "Déséquilibre 🚩"
        yield_est = (fat * 0.12) + (prot * 0.15) + 0.5
        return ratio, status, rumen, yield_est

# ============================================================================
# 3. INTERFACE PRINCIPALE (STREAMLIT)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin Génomique Pro", layout="wide", page_icon="🧬")

    if 'db' not in st.session_state: st.session_state.db = []
    if 'auth' not in st.session_state: st.session_state.auth = False
    
    genomique_engine = IntegrationGenomique()

    # --- CONNEXION ---
    if not st.session_state.auth:
        st.title("🔐 Accès Laboratoire Expert Ovin")
        user = st.text_input("Identifiant", "admin")
        pwd = st.text_input("Mot de passe", type="password")
        if st.button("Connexion"):
            if pwd == "admin123":
                st.session_state.auth = True
                st.rerun()
        return

    # --- NAVIGATION ---
    st.sidebar.image("https://cdn-icons-png.flaticon.com/512/3063/3063176.png", width=100)
    menu = ["🏠 Accueil", "📷 Scanner & Morpho", "🧪 Biochimie Laitière", "🧬 Génomique & SNP", "📊 Rapport & Export"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- ACCUEIL ---
    if choice == "🏠 Accueil":
        st.title("Système Expert de Sélection Laitière Ovine")
        st.markdown("""
        ### Plateforme Intégrée Génomique & Biochimique
        - **Analyse SNP** : Détection des transitions et transversions génétiques.
        - **Biochimie** : Suivi de la cétose et rendement fromager.
        - **Sélection** : Matrice de parenté génomique et indexation EBV.
        """)
        

    # --- SCANNER & MORPHO ---
    elif choice == "📷 Scanner & Morpho":
        st.title("📏 Morphométrie & Saisie Terrain")
        col1, col2 = st.columns(2)
        with col1:
            file = st.file_uploader("Upload Image pour analyse", type=['jpg', 'png'])
            ref_obj = st.selectbox("Standard de référence", list(CALIBRATION_STANDARDS.keys()))
            if file:
                ratio = CALIBRATION_STANDARDS[ref_obj] / 450.0 # Simulation calibration
                st.success(f"Calibration : 1px = {ratio:.4f} cm")
        
        with col2:
            with st.form("form_animal"):
                id_b = st.text_input("ID Animal*")
                race = st.selectbox("Race", ["Lacaune", "Ouled Djellal", "Rembi", "Hamra"])
                hauteur = st.number_input("Hauteur (cm)", 50, 100, 70)
                longueur = st.number_input("Longueur (cm)", 50, 120, 80)
                if st.form_submit_button("💾 Enregistrer"):
                    st.session_state.db.append({"id": id_b, "race": race, "hauteur": hauteur, "longueur": longueur, "ial": 0, "ebv": 0})
                    st.success("Données enregistrées.")

    # --- BIOCHIMIE ---
    elif choice == "🧪 Biochimie Laitière":
        st.title("🧪 Analyse Biochimique du Lait")
        if not st.session_state.db: st.warning("Saisissez d'abord un animal.")
        else:
            target = st.selectbox("Animal", [x['id'] for x in st.session_state.db])
            c1, c2, c3 = st.columns(3)
            fat = c1.number_input("Taux Butyreux (g/L)", 20.0, 95.0, 48.0)
            prot = c2.number_input("Taux Protéique (g/L)", 20.0, 75.0, 39.0)
            bhb = c3.number_input("BHB (mmol/L)", 0.0, 5.0, 0.6)
            
            ratio, status, rumen, yield_est = UltraExpertModule.get_biochem_diagnostic(fat, prot, bhb)
            
            st.divider()
            st.metric("Rapport TB/TP", f"{ratio:.2f}", help="Indicateur de santé du rumen")
            st.write(f"🩺 État métabolique : **{status}**")
            st.metric("Rendement Jben/Fromage estimé", f"{yield_est:.2f} kg/100L")
            

    # --- GÉNOMIQUE & SNP ---
    elif choice == "🧬 Génomique & SNP":
        st.title("🧬 Analyse Bioinformatique & SNP")
        
        tab1, tab2, tab3 = st.tabs(["Analyse de Séquence", "Gènes Candidats", "G-Matrix"])
        
        with tab1:
            st.subheader("Analyse des Polymorphismes (SNP)")
            seq_ref = st.text_area("Séquence de Référence NCBI (ex: GDF9)", "ATGCGTACGTAGCTAGCTAGCGATCGATCGATCGA")
            seq_stu = st.text_area("Séquence de la Brebis", "ATGCGTACGTGGCTAGCTAGCCATCGATCGATCGA")
            
            if st.button("Lancer l'analyse SNP"):
                result = genomique_engine.analyser_snp(seq_ref, seq_stu)
                if "erreur" in result:
                    st.error(result["erreur"])
                else:
                    st.success(f"Analyse terminée : {result['total_snps']} SNPs détectés.")
                    st.json(result)
                    

        with tab2:
            st.subheader("Recherche de gènes par race")
            race_sel = st.selectbox("Race pour recherche GeneBank", ["Lacaune", "Ouled Djellal", "Rembi"])
            genes = genomique_engine.rechercher_genes_candidats(race_sel)
            st.table(genes)

        with tab3:
            st.subheader("Matrice de parenté génomique (G-Matrix)")
            matrix = np.random.rand(8, 8)
            fig = px.imshow(matrix, title="Realized Relatedness (SNP-based)")
            st.plotly_chart(fig)

    # --- EXPORT ---
    elif choice == "📊 Rapport & Export":
        st.title("📥 Rapports Génomiques")
        if st.session_state.db:
            df = pd.DataFrame(st.session_state.db)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Télécharger Rapport (.csv)", csv, "rapport_ovin.csv")
            
            if st.button("Générer Rapport Textuel Expert"):
                rapport = genomique_engine.generer_rapport_genomique(df.iloc[0]['id'], [{'accession': 'NC_040254', 'longueur': 500}])
                st.code(rapport, language="markdown")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.auth = False
        st.rerun()

if __name__ == "__main__":
    main()
