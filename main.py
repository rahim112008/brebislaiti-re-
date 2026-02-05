"""
PROJET : EXPERT OVIN DZ PRO (VERSION INTÉGRALE 2026)
Domaine : Sélection génétique, Génomique, Morphométrie & Gestion Laitière
Auteur : rahim LABORATOIRE GenApAgiE 
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from datetime import datetime
import io
import time
import base64

# ============================================================================
# 1. CONFIGURATION, STANDARDS ET ROLES
# ============================================================================

ROLES = {
    'eleveur': {'label': 'Éleveur', 'perms': ['saisie', 'voir_own', 'excel']},
    'tech': {'label': 'Technicien Conseil', 'perms': ['saisie', 'voir_all', 'stats', 'genetique', 'excel']},
    'admin': {'label': 'Administrateur', 'perms': ['all', 'dashboard', 'users']}
}

CALIBRATION_STANDARDS = {
    "Standard 1m (Mètre ruban)": 100.0,
    "Feuille A4 (Hauteur)": 29.7,
    "Carte Bancaire (Largeur)": 8.56,
    "Pièce 100 DA (Diamètre)": 2.95
}

RACES_DZ = ["Ouled Djellal", "Rembi", "Hamra", "Berbère", "Tergui", "Lacaune"]

# ============================================================================
# 2. MODULE GÉNOMIQUE ET BIO-INFORMATIQUE (GenomicsLab)
# ============================================================================

class GenomicsLab:
    @staticmethod
    def calculate_hwe(p):
        """Équilibre de Hardy-Weinberg : p^2 + 2pq + q^2 = 1"""
        q = 1 - p
        return {"AA": p**2, "Aa": 2*p*q, "aa": q**2}

    @staticmethod
    def verify_filiation(id_a, id_p, id_m):
        """Simulation d'analyse de filiation par microsatellites (STR)"""
        match = np.random.choice([True, False], p=[0.98, 0.02])
        return match

    @staticmethod
    def calculate_ebv(prod_mere, prod_gm, moyenne_pop=250):
        """Calcul simplifié de la Valeur Génétique Estimée (EBV) - Héritabilité h2=0.25"""
        h2 = 0.25
        ebv = h2 * (0.5 * (prod_mere - moyenne_pop) + 0.25 * (prod_gm - moyenne_pop))
        return round(ebv, 2)

# ============================================================================
# 3. MODULE SCANNER IA ET CALIBRATION (AIScanner)
# ============================================================================

class AIScanner:
    @staticmethod
    def process_calibration(image, ref_type):
        """Calcule le ratio pixel/cm basé sur l'objet de référence choisi"""
        ref_reel_cm = CALIBRATION_STANDARDS[ref_type]
        # Simulation d'une détection OpenCV (en production : cv2.findContours)
        pixel_reference_detecte = 450.0 
        ratio = ref_reel_cm / pixel_reference_detecte
        return ratio

# ============================================================================
# 4. MODULE GESTION EXCEL (ExcelManager)
# ============================================================================

class ExcelManager:
    @staticmethod
    def to_excel(df):
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Données_Elevage')
            # Stylisation simplifiée
            workbook = writer.book
            worksheet = writer.sheets['Données_Elevage']
            header_format = workbook.add_format({'bold': True, 'bg_color': '#366092', 'font_color': 'white'})
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
        return output.getvalue()

# ============================================================================
# 5. LOGIQUE D'INTERFACE PRINCIPALE (UI)
# ============================================================================

def main():
    st.set_page_config(page_title="Expert Ovin DZ Pro", layout="wide", page_icon="🐑")

    # Initialisation de la base de données de session
    if 'db' not in st.session_state:
        st.session_state.db = []
    if 'auth' not in st.session_state:
        st.session_state.auth = False

    # --- ÉCRAN DE CONNEXION ---
    if not st.session_state.auth:
        st.title("🔐 Connexion - Expert Ovin DZ Pro")
        col1, col2 = st.columns(2)
        with col1:
            email = st.text_input("Email", "admin@demo.com")
            pwd = st.text_input("Mot de passe", type="password", value="admin123")
            if st.button("Se connecter"):
                if email == "admin@demo.com" and pwd == "admin123":
                    st.session_state.auth = True
                    st.session_state.user_role = "admin"
                    st.session_state.user_name = "Administrateur Central"
                    st.rerun()
                elif email == "eleveur@demo.com" and pwd == "demo123":
                    st.session_state.auth = True
                    st.session_state.user_role = "eleveur"
                    st.session_state.user_name = "Ali Benali (Éleveur)"
                    st.rerun()
                else:
                    st.error("Identifiants incorrects.")
        return

    # --- BARRE LATÉRALE DE NAVIGATION ---
    st.sidebar.title("🐑 Expert Ovin Pro")
    st.sidebar.write(f"Utilisateur : {st.session_state.user_name}")
    menu = ["🏠 Accueil", "📷 Scanner & Saisie", "🧬 Labo Génétique", "⚖️ Comparatif", "📊 Dashboard Admin", "📥 Excel"]
    choice = st.sidebar.radio("Navigation", menu)

    # --- CONTENU DES MODULES ---

    if choice == "🏠 Accueil":
        st.title("Bienvenue sur la plateforme Expert Ovin DZ")
        st.markdown("""
        ### Système Intégré de Sélection Génomique et Morphologique
        Cette application est dédiée à l'amélioration du cheptel Algérien, avec un focus particulier sur les **brebis laitières**.
        
        - **Génétique** : Calcul des index EBV et validation de filiation.
        - **Morphométrie** : Scanner IA calibré (A4, 100DA, Carte).
        - **Phénotype** : Pointage complet de la mamelle.
        """)
        

    elif choice == "📷 Scanner & Saisie":
        st.title("📏 Saisie Morphométrique & Scanner IA")
        
        with st.expander("📷 Option Scanner IA (Calibration)", expanded=True):
            col_img, col_res = st.columns(2)
            with col_img:
                uploaded_file = st.file_uploader("Prendre/Charger photo", type=['jpg', 'jpeg', 'png'])
                ref_type = st.selectbox("Objet de référence pour calibration", list(CALIBRATION_STANDARDS.keys()))
            
            if uploaded_file:
                ratio = AIScanner.process_calibration(uploaded_file, ref_type)
                h_detectee = 700 * ratio # Simulation
                l_detectee = 820 * ratio # Simulation
                with col_res:
                    st.metric("Hauteur (Calibrée)", f"{h_detectee:.1f} cm")
                    st.metric("Longueur (Calibrée)", f"{l_detectee:.1f} cm")
                    st.info(f"Ratio : 1 pixel = {ratio:.4f} cm")

        st.divider()
        
        with st.form("form_saisie"):
            st.subheader("📝 Fiche de l'Animal")
            c1, c2, c3 = st.columns(3)
            id_animal = c1.text_input("ID Animal (Boucle)*")
            race = c2.selectbox("Race", RACES_DZ)
            age = c3.number_input("Âge (mois)", 6, 180, 24)

            st.subheader("📐 Mensurations & Mamelle (Spécial Lait)")
            
            m1, m2, m3 = st.columns(3)
            hauteur = m1.number_input("Hauteur au garrot (cm)", 30.0, 110.0, 70.0)
            longueur = m2.number_input("Longueur de corps (cm)", 40.0, 130.0, 80.0)
            profondeur_mam = m3.slider("Profondeur Mamelle (Score 1-9)", 1, 9, 5)

            st.subheader("🧬 Caractères Phénotypiques")
            p1, p2, p3 = st.columns(3)
            attache = p1.slider("Attache Arrière", 1, 9, 5)
            trayons = p2.slider("Orientation Trayons", 1, 9, 5)
            sillon = p3.slider("Sillon Médian", 1, 9, 5)

            if st.form_submit_button("💾 Enregistrer dans la Base"):
                # Calculs automatiques
                ial = (attache * 0.4) + (trayons * 0.3) + (sillon * 0.3)
                poids = (hauteur**2 * longueur) / 10800
                
                entry = {
                    "id": id_animal, "race": race, "age": age, 
                    "hauteur": hauteur, "longueur": longueur, "poids_est": round(poids, 1),
                    "ial": round(ial, 2), "prof_mamelle": profondeur_mam,
                    "date": datetime.now().strftime("%d/%m/%Y"), "eleveur": st.session_state.user_name
                }
                st.session_state.db.append(entry)
                st.success(f"Animal {id_animal} enregistré avec succès !")

    elif choice == "🧬 Labo Génétique":
        st.title("🔬 Laboratoire de Génétique & Bio-informatique")
        
        tab1, tab2 = st.tabs(["Filiation ADN", "Génétique de Population"])
        
        with tab1:
            st.subheader("🧬 Test de Paternité / Maternité")
            col_f1, col_f2, col_f3 = st.columns(3)
            id_a = col_f1.text_input("ID Agneau")
            id_p = col_f2.text_input("ID Père (Bélier)")
            id_m = col_f3.text_input("ID Mère (Brebis)")
            
            if st.button("Lancer l'analyse de filiation"):
                res = GenomicsLab.verify_filiation(id_a, id_p, id_m)
                if res:
                    st.success(f"✅ Filiation confirmée pour {id_a}")
                else:
                    st.error("❌ Exclusion de paternité détectée.")
                
        with tab2:
            st.subheader("📊 Équilibre de Hardy-Weinberg")
            p = st.slider("Fréquence de l'allèle dominant (p)", 0.0, 1.0, 0.6)
            hwe = GenomicsLab.calculate_hwe(p)
            st.write(f"Distribution génotypique : AA: {hwe['AA']:.2%}, Aa: {hwe['Aa']:.2%}, aa: {hwe['aa']:.2%}")
            
            # Graphique Population
            fig = px.bar(x=['AA', 'Aa', 'aa'], y=[hwe['AA'], hwe['Aa'], hwe['aa']], title="Répartition des Génotypes")
            st.plotly_chart(fig)

    elif choice == "⚖️ Comparatif":
        st.title("⚖️ Mode Comparatif - Sélection de Reproducteurs")
        if len(st.session_state.db) < 2:
            st.warning("Veuillez saisir au moins 2 animaux.")
        else:
            ids = [x['id'] for x in st.session_state.db]
            c1, c2 = st.columns(2)
            with c1:
                sel1 = st.selectbox("Animal A", ids, key="s1")
                data1 = next(i for i in st.session_state.db if i['id'] == sel1)
                st.metric("Index Laitier (IAL)", data1['ial'])
                st.metric("Poids Est.", f"{data1['poids_est']} kg")
            with c2:
                sel2 = st.selectbox("Animal B", ids, key="s2")
                data2 = next(i for i in st.session_state.db if i['id'] == sel2)
                st.metric("Index Laitier (IAL)", data2['ial'], delta=round(data2['ial']-data1['ial'], 2))
                st.metric("Poids Est.", f"{data2['poids_est']} kg", delta=round(data2['poids_est']-data1['poids_est'], 1))

            # Radar Chart
            cat = ['Hauteur', 'Longueur', 'IAL', 'Prof. Mamelle']
            fig = go.Figure()
            fig.add_trace(go.Scatterpolar(r=[data1['hauteur'], data1['longueur'], data1['ial']*10, data1['prof_mamelle']*10], theta=cat, fill='toself', name=sel1))
            fig.add_trace(go.Scatterpolar(r=[data2['hauteur'], data2['longueur'], data2['ial']*10, data2['prof_mamelle']*10], theta=cat, fill='toself', name=sel2))
            st.plotly_chart(fig)

    elif choice == "📊 Dashboard Admin":
        st.title("📊 Tableau de Bord Administrateur")
        if not st.session_state.db:
            st.info("Aucune donnée disponible.")
        else:
            df = pd.DataFrame(st.session_state.db)
            st.dataframe(df, use_container_width=True)
            
            col_d1, col_d2 = st.columns(2)
            with col_d1:
                st.plotly_chart(px.pie(df, names='race', title="Répartition par Race"))
            with col_d2:
                st.plotly_chart(px.box(df, x='race', y='ial', title="Variation de l'Index Laitier"))

    elif choice == "📥 Excel":
        st.title("📥 Import / Export de Données")
        if st.session_state.db:
            df_export = pd.DataFrame(st.session_state.db)
            xls_data = ExcelManager.to_excel(df_export)
            st.download_button(label="📥 Télécharger Rapport Excel Pro", data=xls_data, file_name="export_ovin_dz.xlsx")
        else:
            st.warning("Base vide.")

    if st.sidebar.button("🚪 Déconnexion"):
        st.session_state.auth = False
        st.rerun()

if __name__ == "__main__":
    main()
