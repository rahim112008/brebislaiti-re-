import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
from datetime import datetime, date
from contextlib import contextmanager
import requests  # Pour appeler l'API backend
import json

# ==========================================
# CONFIGURATION ET BASE DE DONNÉES
# ==========================================
DB_NAME = "expert_ovin_simulation.db"

@contextmanager
def get_db_connection():
    conn = sqlite3.connect(DB_NAME, check_same_thread=False)
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()

def init_database():
    with get_db_connection() as conn:
        conn.executescript('''
            CREATE TABLE IF NOT EXISTS animaux (
                id TEXT PRIMARY KEY, boucle TEXT, race TEXT, date_naiss DATE);
            CREATE TABLE IF NOT EXISTS production (
                id INTEGER PRIMARY KEY AUTOINCREMENT, animal_id TEXT, 
                tb REAL, tp REAL, volume_lait REAL,
                FOREIGN KEY(animal_id) REFERENCES animaux(id));
            CREATE TABLE IF NOT EXISTS morphometrie (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                animal_id TEXT NOT NULL,
                date_mesure TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                longueur_corps_cm REAL,
                hauteur_garrot_cm REAL,
                perimetre_thoracique_cm REAL,
                poids_estime_kg REAL,
                FOREIGN KEY (animal_id) REFERENCES animaux(id)
            );
        ''')

# ==========================================
# GÉNÉRATEUR DE POPULATION (20 BREBIS)
# ==========================================
def generer_troupeau_demo(n=20):
    with get_db_connection() as conn:
        cursor = conn.cursor()
        # On vérifie si les données existent déjà
        if cursor.execute("SELECT COUNT(*) FROM animaux").fetchone()[0] == 0:
            for i in range(1, n + 1):
                a_id = f"SIM_{i:02d}"
                boucle = f"FR-161-{i:03d}"
                
                # Simulation de données biologiques (Distribution normale)
                # Moyenne TB: 70g/L, TP: 55g/L, Volume: 1.8L/jour
                tb = round(np.random.normal(70, 5), 2)
                tp = round(np.random.normal(55, 3), 2)
                vol = round(np.random.normal(1.8, 0.4), 2)
                
                cursor.execute("INSERT INTO animaux VALUES (?,?,?,?)", 
                              (a_id, boucle, "Lacaune", "2024-03-15"))
                cursor.execute("INSERT INTO production (animal_id, tb, tp, volume_lait) VALUES (?,?,?,?)", 
                              (a_id, tb, tp, vol))
            return True
    return False

# ==========================================
# MODULE DE VISUALISATION
# ==========================================
def view_production_analysis():
    st.title("📊 Analyse de Production du Troupeau (Simulation n=20)")
    
    with get_db_connection() as conn:
        df = pd.read_sql('''
            SELECT a.boucle, p.tb, p.tp, p.volume_lait 
            FROM animaux a JOIN production p ON a.id = p.animal_id
        ''', conn)

    if not df.empty:
        # Métriques Globales
        c1, c2, c3 = st.columns(3)
        c1.metric("Moyenne TB (Gras)", f"{df['tb'].mean():.2f} g/L")
        c2.metric("Moyenne TP (Protéines)", f"{df['tp'].mean():.2f} g/L")
        c3.metric("Volume Moyen", f"{df['volume_lait'].mean():.2f} L/j")

        # Graphiques Comparatifs
        st.subheader("📈 Comparaison Individuelle : Volume vs Taux Butyreux")
        fig = px.scatter(df, x="volume_lait", y="tb", color="tp",
                         size="tb", hover_name="boucle",
                         title="Analyse Multi-paramètres (Volume, TB, TP)",
                         labels={"volume_lait": "Volume (L)", "tb": "Taux Butyreux (g/L)"})
        st.plotly_chart(fig, use_container_width=True)
        
        st.subheader("📋 Données Détaillées")
        st.dataframe(df.sort_values(by="volume_lait", ascending=False), use_container_width=True)
    else:
        st.info("Cliquez sur le bouton dans la barre latérale pour générer les données.")

# ==========================================
# MODULE SCANNER MORPHO
# ==========================================
def scanner_morpho_module():
    st.title("🐑 Scanner Morphométrique 3D")
    
    # Initialisation de l'état de session pour les photos
    if 'scan_photos' not in st.session_state:
        st.session_state['scan_photos'] = []
    if 'scan_results' not in st.session_state:
        st.session_state['scan_results'] = None
    if 'reference_length' not in st.session_state:
        st.session_state['reference_length'] = 50.0  # Valeur par défaut (cm)
    
    # Étape 1: Capture des photos
    st.header("1️⃣ Capture des photos")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        **Instructions de capture :**
        - Prenez **20-30 photos** du mouton sous différents angles
        - Tournez autour de l'animal de manière régulière
        - Assurez-vous que l'animal est immobile
        - Évitez les reflets et ombres marquées
        """)
        
        # Widget de capture photo
        photo = st.camera_input("Prenez une photo", key="camera_input")
        
        if photo:
            # Ajouter la photo à la liste
            st.session_state['scan_photos'].append(photo.getvalue())
            st.success(f"Photo {len(st.session_state['scan_photos'])} capturée !")
            st.rerun()
    
    with col2:
        st.metric("Photos capturées", f"{len(st.session_state['scan_photos'])}/30")
        
        if st.button("🔄 Réinitialiser", use_container_width=True):
            st.session_state['scan_photos'] = []
            st.session_state['scan_results'] = None
            st.rerun()
    
    # Étape 2: Configuration de l'échelle
    if len(st.session_state['scan_photos']) > 0:
        st.header("2️⃣ Configuration de l'échelle")
        
        st.session_state['reference_length'] = st.number_input(
            "Longueur de référence connue (en cm)",
            min_value=10.0,
            max_value=200.0,
            value=st.session_state['reference_length'],
            help="Mesurez une distance précise sur l'animal (ex: hauteur au garrot)"
        )
        
        st.info(f"🔹 Utilisez un mètre ruban pour mesurer une distance précise sur l'animal.")
    
    # Étape 3: Traitement 3D
    if len(st.session_state['scan_photos']) >= 10:
        st.header("3️⃣ Traitement 3D")
        
        if st.button("🚀 Lancer la reconstruction 3D", type="primary", use_container_width=True):
            with st.spinner("Reconstruction 3D en cours... Cette opération peut prendre 2-3 minutes."):
                try:
                    # Préparer les données pour l'envoi
                    photos_data = st.session_state['scan_photos']
                    reference_length = st.session_state['reference_length']
                    
                    # Simulation pour démonstration
                    results = simulate_3d_processing(reference_length)
                    
                    # Stocker les résultats
                    st.session_state['scan_results'] = results
                    
                    st.success("✅ Analyse 3D terminée avec succès !")
                    
                except Exception as e:
                    st.error(f"❌ Erreur lors du traitement : {str(e)}")
    
    # Étape 4: Affichage des résultats
    if st.session_state['scan_results']:
        st.header("4️⃣ Résultats morphométriques")
        
        results = st.session_state['scan_results']
        
        # Affichage des mesures
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Longueur corps", f"{results['longueur_corps_cm']:.1f} cm")
        with col2:
            st.metric("Hauteur garrot", f"{results['hauteur_garrot_cm']:.1f} cm")
        with col3:
            st.metric("Périmètre thorax", f"{results['perimetre_thoracique_cm']:.1f} cm")
        with col4:
            st.metric("Poids estimé", f"{results['poids_estime_kg']:.1f} kg")
        
        # Graphique de visualisation
        st.subheader("Visualisation des proportions")
        
        metrics_df = pd.DataFrame({
            'Mesure': ['Longueur corps', 'Hauteur garrot', 'Périmètre thorax'],
            'Valeur (cm)': [
                results['longueur_corps_cm'],
                results['hauteur_garrot_cm'],
                results['perimetre_thoracique_cm']
            ]
        })
        
        fig = px.bar(metrics_df, x='Mesure', y='Valeur (cm)', 
                     title="Dimensions corporelles",
                     color='Mesure')
        st.plotly_chart(fig, use_container_width=True)
        
        # Sauvegarde dans la base de données
        st.subheader("💾 Sauvegarde des résultats")
        
        with get_db_connection() as conn:
            animaux_df = pd.read_sql("SELECT id, boucle FROM animaux", conn)
        
        if not animaux_df.empty:
            selected_boucle = st.selectbox(
                "Sélectionnez la boucle de l'animal",
                animaux_df['boucle']
            )
            
            if st.button("Enregistrer les mesures", type="primary"):
                animal_id = animaux_df.loc[animaux_df['boucle'] == selected_boucle, 'id'].iloc[0]
                
                with get_db_connection() as conn:
                    cursor = conn.cursor()
                    cursor.execute('''
                        INSERT INTO morphometrie 
                        (animal_id, longueur_corps_cm, hauteur_garrot_cm, 
                         perimetre_thoracique_cm, poids_estime_kg)
                        VALUES (?, ?, ?, ?, ?)
                    ''', (
                        animal_id,
                        results['longueur_corps_cm'],
                        results['hauteur_garrot_cm'],
                        results['perimetre_thoracique_cm'],
                        results['poids_estime_kg']
                    ))
                
                st.success(f"Mesures enregistrées pour la brebis {selected_boucle} !")
                
                # Afficher l'historique des mesures
                st.subheader("📊 Historique des mesures")
                historique = pd.read_sql(f'''
                    SELECT date_mesure, longueur_corps_cm, hauteur_garrot_cm, 
                           perimetre_thoracique_cm, poids_estime_kg
                    FROM morphometrie 
                    WHERE animal_id = '{animal_id}'
                    ORDER BY date_mesure DESC
                ''', conn)
                
                if not historique.empty:
                    st.dataframe(historique, use_container_width=True)
                else:
                    st.info("Première mesure enregistrée pour cet animal.")
        else:
            st.warning("Générez d'abord le troupeau de démonstration dans la barre latérale.")

# ==========================================
# FONCTIONS DE TRAITEMENT 3D (SIMULATION)
# ==========================================
def simulate_3d_processing(reference_length):
    """Fonction de simulation pour les tests"""
    import random
    
    # Génération de mesures réalistes pour un mouton
    return {
        'longueur_corps_cm': round(reference_length * random.uniform(1.4, 1.8), 1),
        'hauteur_garrot_cm': round(reference_length * random.uniform(1.1, 1.3), 1),
        'perimetre_thoracique_cm': round(reference_length * random.uniform(2.0, 2.4), 1),
        'poids_estime_kg': round((reference_length ** 3) * random.uniform(0.002, 0.003), 1),
        'precision_estimee': round(random.uniform(95, 98), 1)
    }

# ==========================================
# MODULE BIOINFORMATIQUE (EXISTANT)
# ==========================================
def bioinformatique_module():
    st.title("🧬 Module Bioinformatique")
    st.write("Utilisez ce module pour lier la production aux variants génétiques.")
    
    # Votre code existant pour NCBI ici...
    st.info("Intégration NCBI à venir...")

# ==========================================
# MAIN APP
# ==========================================
def main():
    # Configuration de la page
    st.set_page_config(
        page_title="Expert Ovin 3D",
        page_icon="🐑",
        layout="wide"
    )
    
    # Initialisation de la base de données
    init_database()
    
    # Sidebar
    with st.sidebar:
        st.title("🔬 Expert Ovin 3D")
        st.markdown("---")
        
        # Génération des données de démo
        if st.button("🚀 Générer 20 brebis (Démo)", use_container_width=True):
            if generer_troupeau_demo():
                st.success("20 brebis ajoutées !")
                st.rerun()
            else:
                st.warning("Les données existent déjà.")
        
        st.markdown("---")
        
        # Navigation
        menu = st.radio(
            "Navigation",
            ["📊 Dashboard Production", "🐑 Scanner Morpho", "🧬 Bioinformatique", "⚙️ Paramètres"]
        )
    
    # Affichage du module sélectionné
    if menu == "📊 Dashboard Production":
        view_production_analysis()
    elif menu == "🐑 Scanner Morpho":
        scanner_morpho_module()  # AJOUT IMPORTANT : CE BLOC ÉTAIT MANQUANT
    elif menu == "🧬 Bioinformatique":
        bioinformatique_module()
    elif menu == "⚙️ Paramètres":
        st.title("Paramètres")
        st.write("Configuration de l'application...")

# ==========================================
# EXÉCUTION
# ==========================================
if __name__ == "__main__":
    main()
