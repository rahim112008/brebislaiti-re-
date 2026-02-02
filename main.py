import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.express as px
from datetime import datetime, date
from contextlib import contextmanager

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
# MAIN APP
# ==========================================
def main():
    init_database()
    
    st.sidebar.title("🔬 Simulation Expert")
    if st.sidebar.button("🚀 Générer 20 Brebis (Demo)"):
        if generer_troupeau_demo():
            st.sidebar.success("20 brebis ajoutées avec succès !")
        else:
            st.sidebar.warning("Les données existent déjà.")

    menu = st.sidebar.radio("Navigation", ["Dashboard Production", "Bioinformatique", "Scanner Morpho"])

    if menu == "Dashboard Production":
        view_production_analysis()
    elif menu == "Bioinformatique":
        st.title("🧬 Module Bioinfo (NCBI)")
        st.write("Utilisez ce module pour lier la production aux variants génétiques.")
        # Appel de vos fonctions NCBI existantes ici...

if __name__ == "__main__":
    main()
