import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from datetime import datetime, timedelta

# --- CONFIGURATION DE LA PAGE ---
st.set_page_config(page_title="Génétique Ovine DZ - Lab", layout="wide")

# --- MOTEURS DE CALCUL ---
class GeneticLabEngine:
    def calculate_delta_g(self, years, h2=0.30, i=1.2, L=3):
        # ΔG = (i * h² * sigma_p) / L (simplifié)
        gain_annuel = (i * h2 * 0.5) / L
        return [2.0 + (gain_annuel * y) for y in range(years + 1)]

    def align_dna(self, seq1, seq2):
        score = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return (score / len(seq1)) * 100

# lab = GeneticLabEngine()

# --- SIDEBAR (NAVIGATION) ---
st.sidebar.title("🧬 Labo Ovin Élite DZ")
menu = st.sidebar.radio("Navigation", [
    "Tableau de Bord", 
    "Scanner & Morphométrie", 
    "Génomique & Bioinfo", 
    "Biométrie (ACP/ANOVA)", 
    "Reproduction & Éponges",
    "Prédiction Progrès (ΔG)"
])

# --- MODULE 1 : TABLEAU DE BORD ---
if menu == "Tableau de Bord":
    st.header("📊 Centre de Contrôle Génétique")
    col1, col2, col3 = st.columns(3)
    col1.metric("Effectif Total", "450 têtes", "+12%")
    col2.metric("Moyenne Laitière", "2.8 L/j", "+0.4")
    col3.metric("Indice Consanguinité", "4.2%", "-0.5%")
    
    st.subheader("Statut du Troupeau Élite")
    df_data = pd.DataFrame({
        'Race': ['Ouled Djellal', 'Hamra', 'Rumbi', 'Elite DZ'],
        'Performance': [70, 65, 68, 92]
    })
    st.bar_chart(df_data.set_index('Race'))

# --- MODULE 2 : SCANNER & MORPHO (Simulation) ---
elif menu == "Scanner & Morphométrie":
    st.header("📏 Scanner Morphométrique IA")
    st.info("Utilisez l'étalon de 1 mètre pour calibrer les mesures.")
    
    uploaded_file = st.file_uploader("Charger l'image de la brebis", type=['jpg', 'png', 'jpeg'])
    
    col1, col2 = st.columns(2)
    with col1:
        if uploaded_file:
            st.image(uploaded_file, caption="Analyse IA en cours...")
        else:
            st.warning("En attente d'image...")
            
    with col2:
        st.subheader("Saisie des mesures (Calibrage 1m)")
        ht_garrot = st.number_input("Hauteur au garrot (cm)", value=80.0)
        lg_bassin = st.number_input("Largeur Bassin (cm)", value=22.0)
        st.success(f"Ecart au Standard Élite : {ht_garrot - 85} cm")

# --- MODULE 3 : GÉNOMIQUE & BIOINFO ---
elif menu == "Génomique & Bioinfo":
    st.header("🧬 Laboratoire de Bioinformatique")
    
    seq_ref = st.text_area("Séquence ADN Référence (FASTA)", "ATGCGGTACTGA...")
    seq_ind = st.text_area("Séquence Individu Élite", "ATGCGGTACTGT...")
    
    if st.button("Lancer l'Alignement"):
        score = sum(1 for a, b in zip(seq_ref, seq_ind) if a == b) / len(seq_ref) * 100
        st.write(f"**Identité de séquence : {score:.2f}%**")
        st.progress(score / 100)
        
    st.subheader("Détection des SNP (Single Nucleotide Polymorphism)")
    st.table({
        "Marqueur": ["Lact_01", "Heat_Resist", "Milk_Fat"],
        "Génotype": ["A/A", "G/T", "C/C"],
        "Effet": ["Élite", "Modéré", "Supérieur"]
    })

# --- MODULE 4 : BIOMÉTRIE (ACP/ANOVA) ---
elif menu == "Biométrie (ACP/ANOVA)":
    st.header("📈 Analyses Biostatistiques")
    
    st.subheader("Analyse en Composantes Principales (ACP)")
    # Simulation de données pour l'ACP
    pca_data = pd.DataFrame(np.random.randn(50, 2), columns=['Axe Morpho', 'Axe Lait'])
    fig = px.scatter(pca_data, x='Axe Morpho', y='Axe Lait', title="Cartographie Génétique du Troupeau")
    st.plotly_chart(fig)
    
    st.subheader("Héritabilité (h²)")
    st.write("Indice h² calculé pour le caractère 'Lait' : **0.32**")

# --- MODULE 5 : REPRODUCTION ---
elif menu == "Reproduction & Éponges":
    st.header("🐑 Suivi de Reproduction")
    date_pose = st.date_input("Date de pose des éponges", datetime.now())
    
    retrait = date_pose + timedelta(days=14)
    saillie = retrait + timedelta(days=2)
    mise_bas = retrait + timedelta(days=150)
    
    st.warning(f"🔔 Alerte retrait éponge : **{retrait}**")
    st.success(f"📅 Mise-bas prévue : **{mise_bas}**")

# --- MODULE 6 : PROGRÈS GÉNÉTIQUE ---
elif menu == "Prédiction Progrès (ΔG)":
    st.header("🚀 Prédiction du Progrès Génétique")
    years = st.slider("Horizon (Années)", 1, 20, 10)
    h2 = st.slider("Héritabilité (h²)", 0.1, 0.5, 0.3)
    
    # Calcul
    engine = GeneticLabEngine()
    progres = engine.calculate_delta_g(years, h2=h2)
    
    fig_delta = px.line(x=list(range(years+1)), y=progres, 
                        labels={'x': 'Années', 'y': 'Production Lait (L/j)'},
                        title="Evolution de la Race Élite")
    st.plotly_chart(fig_delta)
