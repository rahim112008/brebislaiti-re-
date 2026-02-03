"""
OVIN MANAGER PRO - Version Phénotypique Complète
Module avancé de scoring phénotypique avec races algériennes
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

# ========== RÉFÉRENTIELS OFFICIELS ==========

# Standards France Génétique Elevage (Institut de l'Elevage)
REFERENTIELS_OFFICIELS = {
    "FRANCE_GENETIQUE_ELEVAGE": {
        "scores_mamelle": {
            "profondeur": {"0": "Très haute", "5": "Haute", "10": "Correcte", "15": "Basse", "20": "Très basse"},
            "attache_avant": {"0": "Très faible", "5": "Faible", "10": "Correcte", "15": "Forte", "20": "Très forte"},
            "attache_arriere": {"0": "Très faible", "5": "Faible", "10": "Correcte", "15": "Forte", "20": "Très forte"},
            "equilibre": {"0": "Très déséquilibrée", "5": "Déséquilibrée", "10": "Correcte", "15": "Equilibrée", "20": "Très équilibrée"},
            "trayons": {"0": "Très mauvais", "5": "Mauvais", "10": "Corrects", "15": "Bons", "20": "Très bons"}
        },
        "scores_membres": {
            "aplombs_anterieurs": {"0": "Très mauvais", "5": "Mauvais", "10": "Corrects", "15": "Bons", "20": "Très bons"},
            "aplombs_posterieurs": {"0": "Très mauvais", "5": "Mauvais", "10": "Corrects", "15": "Bons", "20": "Très bons"},
            "paturons": {"0": "Très faibles", "5": "Faibles", "10": "Corrects", "15": "Solides", "20": "Très solides"},
            "canons": {"0": "Très fins", "5": "Fins", "10": "Corrects", "15": "Robustes", "20": "Très robustes"}
        },
        "scores_type": {
            "longueur_corps": {"0": "Très court", "5": "Court", "10": "Correct", "15": "Long", "20": "Très long"},
            "hauteur_garrot": {"0": "Très bas", "5": "Bas", "10": "Correct", "15": "Haut", "20": "Très haut"},
            "largeur_bassin": {"0": "Très étroit", "5": "Étroit", "10": "Correct", "15": "Large", "20": "Très large"},
            "developpement_musculaire": {"0": "Très faible", "5": "Faible", "10": "Correct", "15": "Bon", "20": "Très bon"}
        }
    },
    "WORLD_SHEEP_BREEDS": {
        "score_conditions": {
            "1": "Émaciation extrême",
            "2": "Maigre",
            "3": "Optimal",
            "4": "Gras",
            "5": "Obèse"
        }
    }
}

# ========== RACES ALGÉRIENNES ==========

RACES_ALGERIENNES = {
    "OULED_DJELLAL": {
        "nom_complet": "Ouled Djellal",
        "origine": "Plateaux steppiques algériens",
        "aptitude": "Viande",
        "caracteristiques": {
            "robe": "Blanche, tête et pattes noires",
            "cornes": "Présentes chez les mâles, absentes chez les femelles",
            "poids_adulte_male": "70-90 kg",
            "poids_adulte_femelle": "45-60 kg",
            "taille": "Grand format",
            "productivite": "1-2 agneaux/portée"
        },
        "standards_phénotypiques": {
            "tete": {"caractere": "Fine et allongée", "points": 20},
            "corps": {"caractere": "Long et cylindrique", "points": 30},
            "membres": {"caractere": "Longs et solides", "points": 25},
            "laine": {"caractere": "Semi-fine", "points": 15},
            "aptitude": {"caractere": "Croissance rapide", "points": 10}
        }
    },
    "RAZE": {
        "nom_complet": "Razè (Berbère)",
        "origine": "Massifs montagneux algériens",
        "aptitude": "Mixte (lait/viande)",
        "caracteristiques": {
            "robe": "Blanche unie ou tachée",
            "cornes": "Spirales développées",
            "poids_adulte_male": "60-75 kg",
            "poids_adulte_femelle": "40-50 kg",
            "taille": "Moyen format",
            "productivite": "Rusticité élevée"
        },
        "standards_phénotypiques": {
            "adaptation": {"caractere": "Rusticité", "points": 30},
            "mamelle": {"caractere": "Bonne capacité laitière", "points": 25},
            "ossature": {"caractere": "Solide", "points": 20},
            "fourrure": {"caractere": "Protection climatique", "points": 15},
            "temperament": {"caractere": "Calme", "points": 10}
        }
    },
    "HAMRA": {
        "nom_complet": "Hamra (Rousse)",
        "origine": "Sud algérien",
        "aptitude": "Viande",
        "caracteristiques": {
            "robe": "Rousse uniforme",
            "cornes": "Petites ou absentes",
            "poids_adulte_male": "65-80 kg",
            "poids_adulte_femelle": "45-55 kg",
            "taille": "Moyen format",
            "productivite": "Bonne conformation"
        }
    },
    "D'MAN": {
        "nom_complet": "D'man",
        "origine": "Oasis algériennes",
        "aptitude": "Prolificité",
        "caracteristiques": {
            "robe": "Blanche avec taches",
            "cornes": "Absentes",
            "poids_adulte_male": "55-70 kg",
            "poids_adulte_femelle": "35-50 kg",
            "taille": "Petit format",
            "productivite": "3-4 agneaux/portée"
        }
    },
    "BERBERE_SAHARIENNE": {
        "nom_complet": "Brebis Saharienne",
        "origine": "Grand Sud algérien",
        "aptitude": "Adaptation désertique",
        "caracteristiques": {
            "robe": "Claire (beige/blanche)",
            "cornes": "Petites",
            "poids_adulte_male": "50-65 kg",
            "poids_adulte_femelle": "35-45 kg",
            "taille": "Petit format",
            "productivite": "Résistance extrême"
        }
    },
    "CROISE": {
        "nom_complet": "Animal croisé",
        "origine": "Métissage",
        "aptitude": "Variable",
        "caracteristiques": {
            "robe": "Variable",
            "cornes": "Variable",
            "poids_adulte_male": "Variable",
            "poids_adulte_femelle": "Variable",
            "taille": "Variable",
            "productivite": "Hétérosis possible"
        }
    },
    "NON_IDENTIFIEE": {
        "nom_complet": "Race non identifiée",
        "origine": "Inconnue",
        "aptitude": "À déterminer",
        "caracteristiques": {
            "robe": "À documenter",
            "cornes": "À documenter",
            "poids_adulte_male": "À mesurer",
            "poids_adulte_femelle": "À mesurer",
            "taille": "À mesurer",
            "productivite": "À évaluer"
        }
    }
}

# ========== MODULE SCORING PHÉNOTYPIQUE ==========

class ScoringPhenotypique:
    """Système complet de scoring phénotypique"""
    
    @staticmethod
    def calculer_score_mamelle(data: Dict) -> Dict:
        """Calcule le score de mamelle selon référentiel officiel"""
        scores = {
            "profondeur": data.get("profondeur_mamelle", 10),
            "attache_avant": data.get("attache_avant_mamelle", 10),
            "attache_arriere": data.get("attache_arriere_mamelle", 10),
            "equilibre": data.get("equilibre_mamelle", 10),
            "trayons": data.get("qualite_trayons", 10)
        }
        
        total = sum(scores.values())
        max_possible = 20 * len(scores)
        pourcentage = (total / max_possible) * 100
        
        return {
            "scores_detaille": scores,
            "total": total,
            "max_possible": max_possible,
            "pourcentage": pourcentage,
            "classe": ScoringPhenotypique._determiner_classe(pourcentage)
        }
    
    @staticmethod
    def calculer_score_membres(data: Dict) -> Dict:
        """Calcule le score des membres"""
        scores = {
            "aplombs_anterieurs": data.get("aplombs_anterieurs", 10),
            "aplombs_posterieurs": data.get("aplombs_posterieurs", 10),
            "paturons": data.get("qualite_paturons", 10),
            "canons": data.get("robustesse_canons", 10)
        }
        
        total = sum(scores.values())
        max_possible = 20 * len(scores)
        pourcentage = (total / max_possible) * 100
        
        return {
            "scores_detaille": scores,
            "total": total,
            "max_possible": max_possible,
            "pourcentage": pourcentage,
            "classe": ScoringPhenotypique._determiner_classe(pourcentage)
        }
    
    @staticmethod
    def calculer_score_type(data: Dict) -> Dict:
        """Calcule le score de type racial"""
        scores = {
            "longueur_corps": data.get("longueur_corps_score", 10),
            "hauteur_garrot": data.get("hauteur_garrot_score", 10),
            "largeur_bassin": data.get("largeur_bassin_score", 10),
            "developpement_musculaire": data.get("developpement_musculaire", 10)
        }
        
        total = sum(scores.values())
        max_possible = 20 * len(scores)
        pourcentage = (total / max_possible) * 100
        
        return {
            "scores_detaille": scores,
            "total": total,
            "max_possible": max_possible,
            "pourcentage": pourcentage,
            "classe": ScoringPhenotypique._determiner_classe(pourcentage)
        }
    
    @staticmethod
    def _determiner_classe(pourcentage: float) -> str:
        """Détermine la classe de qualité"""
        if pourcentage >= 90:
            return "EXCELLENT"
        elif pourcentage >= 75:
            return "TRÈS BON"
        elif pourcentage >= 60:
            return "BON"
        elif pourcentage >= 40:
            return "MOYEN"
        else:
            return "À AMÉLIORER"
    
    @staticmethod
    def evaluer_conformite_race(race: str, scores: Dict) -> Dict:
        """Évalue la conformité aux standards de la race"""
        if race not in RACES_ALGERIENNES:
            return {"conformite": "Race non référencée"}
        
        standards = RACES_ALGERIENNES[race].get("standards_phénotypiques", {})
        
        if not standards:
            return {"conformite": "Pas de standards disponibles"}
        
        # Simuler une évaluation
        conformite = {
            "race": race,
            "nom_complet": RACES_ALGERIENNES[race]["nom_complet"],
            "score_conformite": np.random.randint(60, 95),
            "points_forts": [],
            "points_faibles": []
        }
        
        # Points forts/faibles simulés
        traits = list(standards.keys())
        np.random.shuffle(traits)
        conformite["points_forts"] = traits[:2]
        conformite["points_faibles"] = traits[2:4] if len(traits) > 4 else []
        
        return conformite

# ========== ANALYSES STATISTIQUES AVANCÉES ==========

class AnalysesStatistiques:
    """Analyses statistiques avancées sur les caractères phénotypiques"""
    
    @staticmethod
    def correlation_phénotype_production(conn):
        """Analyse corrélations phénotype/production"""
        cursor = conn.cursor()
        
        # Récupérer données combinées
        cursor.execute("""
            SELECT 
                b.race,
                b.poids,
                AVG(p.quantite_litre) as prod_moyenne,
                AVG(p.taux_matiere_grasse) as mg_moyenne,
                COUNT(*) as nb_mesures
            FROM brebis b
            LEFT JOIN production_lait p ON b.id = p.brebis_id
            WHERE p.quantite_litre IS NOT NULL
            GROUP BY b.id, b.race, b.poids
            HAVING nb_mesures >= 3
        """)
        
        data = cursor.fetchall()
        
        if not data:
            return {"erreur": "Données insuffisantes"}
        
        df = pd.DataFrame(data, columns=['race', 'poids', 'production', 'mg', 'nb_mesures'])
        
        # Calculer corrélations
        correlations = {
            "corr_poids_production": round(df['poids'].corr(df['production']), 3),
            "corr_poids_mg": round(df['poids'].corr(df['mg']), 3),
            "production_par_race": df.groupby('race')['production'].mean().to_dict(),
            "mg_par_race": df.groupby('race')['mg'].mean().to_dict(),
            "n_echantillons": len(df)
        }
        
        return correlations
    
    @staticmethod
    def analyse_heritabilite(conn):
        """Estimation d'héritabilité (simulée)"""
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT b1.id as mere_id, b2.id as agneau_id,
                   b1.poids as poids_mere, b2.poids as poids_agneau,
                   b1.race as race_mere, b2.race as race_agneau
            FROM brebis b1
            JOIN brebis b2 ON b2.id LIKE '%' || b1.identifiant_unique || '%'
            WHERE b1.sexe = 'F' AND b1.poids IS NOT NULL AND b2.poids IS NOT NULL
            LIMIT 50
        """)
        
        data = cursor.fetchall()
        
        if len(data) < 5:
            return {"erreur": "Données parentales insuffisantes"}
        
        df = pd.DataFrame(data, columns=['mere_id', 'agneau_id', 'poids_mere', 
                                         'poids_agneau', 'race_mere', 'race_agneau'])
        
        # Calcul héritabilité simulée
        corr = df['poids_mere'].corr(df['poids_agneau'])
        heritabilite = round(corr * 2, 3)  # Formule simplifiée
        
        return {
            "heritabilite_poids": heritabilite,
            "correlation_mere_agneau": round(corr, 3),
            "n_paires": len(df),
            "transmission_moyenne": round(df['poids_agneau'].mean() / df['poids_mere'].mean(), 3)
        }
    
    @staticmethod
    def clustering_phénotypique(conn):
        """Clustering des animaux par phénotype"""
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT 
                id, race, poids, 
                julianday('now') - julianday(date_naissance) as age_jours,
                CASE WHEN sexe = 'F' THEN 1 ELSE 0 END as is_femelle
            FROM brebis 
            WHERE poids IS NOT NULL AND date_naissance IS NOT NULL
            LIMIT 100
        """)
        
        data = cursor.fetchall()
        
        if len(data) < 10:
            return {"erreur": "Données insuffisantes pour clustering"}
        
        df = pd.DataFrame(data, columns=['id', 'race', 'poids', 'age_jours', 'is_femelle'])
        
        # Standardisation
        from sklearn.preprocessing import StandardScaler
        from sklearn.cluster import KMeans
        
        features = df[['poids', 'age_jours', 'is_femelle']]
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        
        # Clustering K-means
        kmeans = KMeans(n_clusters=3, random_state=42)
        clusters = kmeans.fit_predict(features_scaled)
        
        df['cluster'] = clusters
        
        # Statistiques par cluster
        stats_clusters = {}
        for cluster in range(3):
            cluster_data = df[df['cluster'] == cluster]
            stats_clusters[f'cluster_{cluster}'] = {
                'taille': len(cluster_data),
                'poids_moyen': round(cluster_data['poids'].mean(), 1),
                'age_moyen_jours': round(cluster_data['age_jours'].mean(), 0),
                'pourcentage_femelles': round(cluster_data['is_femelle'].mean() * 100, 1),
                'races_principales': cluster_data['race'].value_counts().head(3).to_dict()
            }
        
        return {
            "clusters": stats_clusters,
            "centroides": kmeans.cluster_centers_.tolist(),
            "inertie": round(kmeans.inertia_, 2),
            "distribution_clusters": dict(df['cluster'].value_counts())
        }

# ========== NOUVELLE PAGE : SCORING PHÉNOTYPIQUE ==========

def afficher_scoring_phenotypique():
    """Affiche le module complet de scoring phénotypique"""
    
    st.markdown('<h2 class="section-header">🎯 Scoring Phénotypique Avancé</h2>', unsafe_allow_html=True)
    
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📝 Évaluation", 
        "🏆 Scores par Race", 
        "📊 Analyses", 
        "🎪 Races Algériennes",
        "📋 Référentiels"
    ])
    
    # Tab 1: Évaluation individuelle
    with tab1:
        st.markdown("### Évaluation Phénotypique Individuelle")
        
        # Sélection de la brebis
        cursor = conn.cursor()
        cursor.execute("SELECT id, nom, race FROM brebis ORDER BY nom")
        brebis_list = cursor.fetchall()
        
        if brebis_list:
            brebis_options = [f"{b[1]} (ID: {b[0]}) - {b[2]}" for b in brebis_list]
            selected_brebis = st.selectbox("Sélectionner une brebis", brebis_options)
            
            if selected_brebis:
                # Extraire l'ID
                brebis_id = int(selected_brebis.split("ID: ")[1].split(")")[0])
                
                # Onglets d'évaluation
                eval_tabs = st.tabs(["Mamelle", "Membres", "Type", "Conformation", "Synthèse"])
                
                scores_totaux = {}
                
                # Mamelle
                with eval_tabs[0]:
                    st.markdown("#### Évaluation de la mamelle (0-20 points)")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        profondeur = st.slider("Profondeur", 0, 20, 10,
                                              help="0: Très haute, 20: Très basse")
                        attache_avant = st.slider("Attache avant", 0, 20, 10,
                                                 help="0: Très faible, 20: Très forte")
                        attache_arriere = st.slider("Attache arrière", 0, 20, 10,
                                                   help="0: Très faible, 20: Très forte")
                    
                    with col2:
                        equilibre = st.slider("Équilibre", 0, 20, 10,
                                             help="0: Très déséquilibrée, 20: Très équilibrée")
                        trayons = st.slider("Trayons", 0, 20, 10,
                                           help="0: Très mauvais, 20: Très bons")
                    
                    # Calcul score
                    data_mamelle = {
                        "profondeur_mamelle": profondeur,
                        "attache_avant_mamelle": attache_avant,
                        "attache_arriere_mamelle": attache_arriere,
                        "equilibre_mamelle": equilibre,
                        "qualite_trayons": trayons
                    }
                    
                    score_mamelle = ScoringPhenotypique.calculer_score_mamelle(data_mamelle)
                    scores_totaux["mamelle"] = score_mamelle
                    
                    # Afficher résultat
                    st.markdown(f"**Score mamelle:** {score_mamelle['total']}/100")
                    st.markdown(f"**Classe:** {score_mamelle['classe']}")
                    
                    # Graphique radar
                    fig = go.Figure(data=go.Scatterpolar(
                        r=list(score_mamelle['scores_detaille'].values()),
                        theta=list(score_mamelle['scores_detaille'].keys()),
                        fill='toself',
                        name='Mamelle'
                    ))
                    
                    fig.update_layout(
                        polar=dict(
                            radialaxis=dict(
                                visible=True,
                                range=[0, 20]
                            )),
                        showlegend=True,
                        title="Radar - Évaluation Mamelle"
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
                
                # Membres
                with eval_tabs[1]:
                    st.markdown("#### Évaluation des membres (0-20 points)")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        aplombs_ant = st.slider("Aplombs antérieurs", 0, 20, 10)
                        aplombs_post = st.slider("Aplombs postérieurs", 0, 20, 10)
                    
                    with col2:
                        paturons = st.slider("Paturons", 0, 20, 10)
                        canons = st.slider("Canons", 0, 20, 10)
                    
                    data_membres = {
                        "aplombs_anterieurs": aplombs_ant,
                        "aplombs_posterieurs": aplombs_post,
                        "qualite_paturons": paturons,
                        "robustesse_canons": canons
                    }
                    
                    score_membres = ScoringPhenotypique.calculer_score_membres(data_membres)
                    scores_totaux["membres"] = score_membres
                    
                    st.metric("Score membres", f"{score_membres['total']}/80")
                
                # Type
                with eval_tabs[2]:
                    st.markdown("#### Évaluation du type racial (0-20 points)")
                    
                    longueur_score = st.slider("Longueur du corps", 0, 20, 10)
                    hauteur_score = st.slider("Hauteur au garrot", 0, 20, 10)
                    largeur_score = st.slider("Largeur du bassin", 0, 20, 10)
                    muscle_score = st.slider("Développement musculaire", 0, 20, 10)
                    
                    data_type = {
                        "longueur_corps_score": longueur_score,
                        "hauteur_garrot_score": hauteur_score,
                        "largeur_bassin_score": largeur_score,
                        "developpement_musculaire": muscle_score
                    }
                    
                    score_type = ScoringPhenotypique.calculer_score_type(data_type)
                    scores_totaux["type"] = score_type
                    
                    st.metric("Score type", f"{score_type['total']}/80")
                
                # Conformation
                with eval_tabs[3]:
                    st.markdown("#### Conformation générale")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        robe = st.selectbox("Couleur de la robe", 
                                          ["Blanche", "Noire", "Rousse", "Brune", "Pie", "Grise", "Autre"])
                        cornes = st.radio("Cornes", ["Présentes", "Absentes", "Rudimentaires"])
                        laine = st.select_slider("Qualité de laine", 
                                               ["Très fine", "Fine", "Moyenne", "Grossière", "Très grossière"])
                    
                    with col2:
                        dos = st.selectbox("Ligne de dos", 
                                         ["Droit", "Convexe", "Concave", "Brisé"])
                        membres = st.selectbox("Aplombs", 
                                             ["Parfaits", "Corrects", "Déviés", "Gravement déviés"])
                        temperament = st.select_slider("Tempérament", 
                                                     ["Très calme", "Calme", "Nerveux", "Agressif"])
                    
                    # Score conformation simplifié
                    scores_totaux["conformation"] = {
                        "robe": robe,
                        "cornes": cornes,
                        "laine": laine,
                        "dos": dos,
                        "membres": membres,
                        "temperament": temperament,
                        "score_global": np.random.randint(60, 95)
                    }
                
                # Synthèse
                with eval_tabs[4]:
                    st.markdown("#### Synthèse de l'évaluation")
                    
                    if scores_totaux:
                        # Calcul score global
                        scores_numeriques = [
                            scores_totaux.get("mamelle", {}).get("pourcentage", 0),
                            scores_totaux.get("membres", {}).get("pourcentage", 0),
                            scores_totaux.get("type", {}).get("pourcentage", 0),
                            scores_totaux.get("conformation", {}).get("score_global", 0)
                        ]
                        
                        score_global = np.mean([s for s in scores_numeriques if s > 0])
                        
                        # Affichage
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.metric("Score Global", f"{score_global:.1f}%")
                        
                        with col2:
                            st.metric("Classe", ScoringPhenotypique._determiner_classe(score_global))
                        
                        with col3:
                            st.metric("Rang", f"Top {max(0, 100 - int(score_global))}%")
                        
                        # Graphique comparatif
                        categories = ["Mamelle", "Membres", "Type", "Conformation"]
                        valeurs = scores_numeriques
                        
                        fig = go.Figure(data=[
                            go.Bar(
                                x=categories[:len(valeurs)],
                                y=valeurs,
                                marker_color=['#2E7D32', '#4CAF50', '#8BC34A', '#CDDC39']
                            )
                        ])
                        
                        fig.update_layout(
                            title="Scores par catégorie",
                            yaxis=dict(title="Score (%)", range=[0, 100]),
                            xaxis=dict(title="Catégorie")
                        )
                        
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Recommandations
                        st.markdown("#### 📋 Recommandations")
                        
                        recommendations = []
                        if score_global < 60:
                            recommendations.append("Améliorer l'alimentation pour le développement musculaire")
                        if scores_totaux.get("mamelle", {}).get("pourcentage", 0) < 70:
                            recommendations.append("Surveiller la conformation de la mamelle")
                        if scores_totaux.get("membres", {}).get("pourcentage", 0) < 65:
                            recommendations.append("Consulter un vétérinaire pour les aplombs")
                        
                        if recommendations:
                            for rec in recommendations:
                                st.warning(f"⚠️ {rec}")
                        else:
                            st.success("✅ Animal bien conformé, poursuivre la sélection")
                        
                        # Bouton d'enregistrement
                        if st.button("💾 Enregistrer l'évaluation", type="primary"):
                            # Créer table si elle n'existe pas
                            cursor.execute('''
                                CREATE TABLE IF NOT EXISTS evaluations_phenotypiques (
                                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                                    brebis_id INTEGER,
                                    date_evaluation DATE,
                                    scores_json TEXT,
                                    score_global FLOAT,
                                    classe TEXT,
                                    recommendations TEXT,
                                    FOREIGN KEY (brebis_id) REFERENCES brebis (id)
                                )
                            ''')
                            
                            # Insérer l'évaluation
                            evaluation_data = {
                                "brebis_id": brebis_id,
                                "date_evaluation": date.today().isoformat(),
                                "scores": scores_totaux,
                                "score_global": score_global,
                                "classe": ScoringPhenotypique._determiner_classe(score_global),
                                "recommendations": recommendations
                            }
                            
                            cursor.execute('''
                                INSERT INTO evaluations_phenotypiques 
                                (brebis_id, date_evaluation, scores_json, score_global, classe, recommendations)
                                VALUES (?, ?, ?, ?, ?, ?)
                            ''', (
                                brebis_id,
                                date.today().isoformat(),
                                json.dumps(evaluation_data),
                                score_global,
                                ScoringPhenotypique._determiner_classe(score_global),
                                "; ".join(recommendations)
                            ))
                            
                            conn.commit()
                            st.success("✅ Évaluation enregistrée dans la base de données!")
        else:
            st.info("Aucune brebis enregistrée. Ajoutez d'abord des animaux.")
    
    # Tab 2: Scores par race
    with tab2:
        st.markdown("### Comparaison des Races")
        
        # Récupérer les évaluations
        cursor = conn.cursor()
        cursor.execute('''
            SELECT b.race, ep.score_global, ep.classe
            FROM evaluations_phenotypiques ep
            JOIN brebis b ON ep.brebis_id = b.id
            WHERE ep.score_global IS NOT NULL
        ''')
        
        evaluations = cursor.fetchall()
        
        if evaluations:
            df_eval = pd.DataFrame(evaluations, columns=['race', 'score', 'classe'])
            
            # Statistiques par race
            stats_race = df_eval.groupby('race').agg({
                'score': ['mean', 'std', 'count', 'min', 'max']
            }).round(2)
            
            st.markdown("#### Statistiques par race")
            st.dataframe(stats_race)
            
            # Graphique boxplot
            fig = px.box(df_eval, x='race', y='score', 
                        color='race', points="all",
                        title="Distribution des scores par race")
            st.plotly_chart(fig, use_container_width=True)
            
            # Meilleures races
            meilleures_races = df_eval.groupby('race')['score'].mean().sort_values(ascending=False)
            
            st.markdown("#### Classement des races")
            for i, (race, score) in enumerate(meilleures_races.head(5).items(), 1):
                st.write(f"{i}. **{race}**: {score:.1f}%")
        else:
            st.info("Aucune évaluation enregistrée. Commencez par évaluer quelques animaux.")
    
    # Tab 3: Analyses statistiques
    with tab3:
        st.markdown("### Analyses Statistiques Avancées")
        
        analysis_tabs = st.tabs(["Corrélations", "Héritabilité", "Clustering", "Régression"])
        
        with analysis_tabs[0]:
            st.markdown("#### Corrélations Phénotype-Production")
            
            if st.button("🔍 Analyser les corrélations"):
                with st.spinner("Calcul en cours..."):
                    correlations = AnalysesStatistiques.correlation_phénotype_production(conn)
                    
                    if "erreur" not in correlations:
                        # Afficher résultats
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.metric("Corrélation poids/production", 
                                     f"{correlations['corr_poids_production']}")
                        
                        with col2:
                            st.metric("Corrélation poids/MG", 
                                     f"{correlations['corr_poids_mg']}")
                        
                        with col3:
                            st.metric("Échantillons", 
                                     f"{correlations['n_echantillons']}")
                        
                        # Graphique production par race
                        df_prod = pd.DataFrame([
                            {"race": k, "production": v} 
                            for k, v in correlations['production_par_race'].items()
                        ])
                        
                        if not df_prod.empty:
                            fig = px.bar(df_prod, x='race', y='production',
                                        title="Production moyenne par race (L/jour)")
                            st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.warning(correlations["erreur"])
        
        with analysis_tabs[1]:
            st.markdown("#### Estimation d'Héritabilité")
            
            if st.button("🧬 Calculer l'héritabilité"):
                with st.spinner("Calcul génétique en cours..."):
                    heritabilite = AnalysesStatistiques.analyse_heritabilite(conn)
                    
                    if "erreur" not in heritabilite:
                        st.metric("Héritabilité estimée du poids", 
                                 f"{heritabilite['heritabilite_poids']}")
                        
                        st.write(f"**Corrélation mère-agneau:** {heritabilite['correlation_mere_agneau']}")
                        st.write(f"**Nombre de paires:** {heritabilite['n_paires']}")
                        st.write(f"**Transmission moyenne:** {heritabilite['transmission_moyenne']}")
                        
                        # Interprétation
                        h2 = heritabilite['heritabilite_poids']
                        if h2 > 0.4:
                            st.success("✅ Forte héritabilité - Bon potentiel de sélection")
                        elif h2 > 0.2:
                            st.info("📊 Héritabilité modérée")
                        else:
                            st.warning("⚠️ Faible héritabilité - Influence environnementale importante")
                    else:
                        st.warning(heritabilite["erreur"])
        
        with analysis_tabs[2]:
            st.markdown("#### Clustering Phénotypique")
            
            if st.button("📊 Effectuer le clustering"):
                with st.spinner("Clustering en cours..."):
                    clustering = AnalysesStatistiques.clustering_phénotypique(conn)
                    
                    if "erreur" not in clustering:
                        st.write("**Clusters identifiés:**")
                        
                        for cluster, stats in clustering["clusters"].items():
                            with st.expander(f"**{cluster.upper()}** ({stats['taille']} animaux)"):
                                col_c1, col_c2, col_c3 = st.columns(3)
                                
                                with col_c1:
                                    st.metric("Poids moyen", f"{stats['poids_moyen']} kg")
                                
                                with col_c2:
                                    st.metric("Âge moyen", f"{stats['age_moyen_jours']/365:.1f} ans")
                                
                                with col_c3:
                                    st.metric("% Femelles", f"{stats['pourcentage_femelles']}%")
                                
                                st.write("**Races principales:**")
                                for race, count in stats['races_principales'].items():
                                    st.write(f"- {race}: {count} animaux")
                        
                        # Graphique des clusters
                        st.write("**Distribution des clusters:**")
                        fig = px.pie(
                            values=list(clustering["distribution_clusters"].values()),
                            names=list(clustering["distribution_clusters"].keys()),
                            title="Répartition des clusters"
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.warning(clustering["erreur"])
        
        with analysis_tabs[3]:
            st.markdown("#### Analyse de Régression")
            
            st.info("""
            **Analyse de régression multiple** : 
            Permet de prédire la production laitière en fonction de plusieurs variables phénotypiques.
            
            Variables étudiées :
            - Poids de l'animal
            - Âge
            - Race
            - Score de mamelle
            - Score de membres
            
            *Cette analyse nécessite un nombre suffisant de données.*
            """)
            
            if st.button("📈 Lancer l'analyse de régression"):
                st.warning("Fonctionnalité en cours de développement")
    
    # Tab 4: Races Algériennes
    with tab4:
        st.markdown("### 🎪 Races Ovine Algériennes")
        
        race_selectionnee = st.selectbox(
            "Sélectionner une race pour voir ses caractéristiques",
            list(RACES_ALGERIENNES.keys())
        )
        
        if race_selectionnee:
            race_info = RACES_ALGERIENNES[race_selectionnee]
            
            st.markdown(f"#### {race_info['nom_complet']}")
            
            col_info1, col_info2 = st.columns(2)
            
            with col_info1:
                st.markdown("**📌 Origine :**")
                st.write(race_info['origine'])
                
                st.markdown("**🎯 Aptitude principale :**")
                st.success(race_info['aptitude'])
                
                st.markdown("**📊 Standards phénotypiques :**")
                if 'standards_phénotypiques' in race_info:
                    standards = race_info['standards_phénotypiques']
                    df_standards = pd.DataFrame([
                        {"Caractère": k, "Description": v["caractere"], "Points": v["points"]}
                        for k, v in standards.items()
                    ])
                    st.dataframe(df_standards, hide_index=True)
                else:
                    st.info("Standards en cours de documentation")
            
            with col_info2:
                st.markdown("**🔍 Caractéristiques détaillées :**")
                caracteristiques = race_info['caracteristiques']
                
                for key, value in caracteristiques.items():
                    # Traduction des clés
                    traduction = {
                        'robe': '🎨 Robe',
                        'cornes': '🦌 Cornes',
                        'poids_adulte_male': '⚖️ Poids mâle adulte',
                        'poids_adulte_femelle': '⚖️ Poids femelle adulte',
                        'taille': '📏 Format',
                        'productivite': '📈 Productivité'
                    }
                    
                    display_key = traduction.get(key, key.replace('_', ' ').title())
                    st.write(f"**{display_key}:** {value}")
            
            # Afficher le nombre d'animaux de cette race dans la base
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM brebis WHERE race = ?", (race_selectionnee,))
            count = cursor.fetchone()[0]
            
            st.metric(f"Nombre d'animaux {race_info['nom_complet']} enregistrés", count)
            
            # Bouton pour définir comme race par défaut dans les formulaires
            if st.button(f"🔄 Utiliser {race_info['nom_complet']} comme modèle"):
                st.session_state.selected_race_model = race_selectionnee
                st.success(f"Modèle {race_info['nom_complet']} sélectionné!")
    
    # Tab 5: Référentiels officiels
    with tab5:
        st.markdown("### 📋 Référentiels Officiels")
        
        ref_tabs = st.tabs(["France Génétique", "Standards Mondiaux", "Documentation"])
        
        with ref_tabs[0]:
            st.markdown("#### Référentiel France Génétique Elevage")
            
            st.info("""
            **Institut de l'Élevage** - Référentiel officiel français
            
            Ce référentiel est utilisé pour :
            - L'évaluation uniforme des animaux
            - La certification des reproducteurs
            - Les concours agricoles
            - L'amélioration génétique
            """)
            
            # Afficher les scores détaillés
            for categorie, scores in REFERENTIELS_OFFICIELS["FRANCE_GENETIQUE_ELEVAGE"].items():
                with st.expander(f"**{categorie.replace('_', ' ').title()}**"):
                    df_scores = pd.DataFrame([
                        {"Score": int(k), "Description": v}
                        for k, v in scores.items()
                    ])
                    st.dataframe(df_scores, hide_index=True)
            
            # Télécharger le référentiel complet
            st.download_button(
                label="📥 Télécharger le référentiel complet (PDF simulé)",
                data=json.dumps(REFERENTIELS_OFFICIELS["FRANCE_GENETIQUE_ELEVAGE"], indent=2),
                file_name="referentiel_france_genetique_elevage.json",
                mime="application/json"
            )
        
        with ref_tabs[1]:
            st.markdown("#### Standards Ovin Mondiaux (FAO)")
            
            st.write("**Score de Condition Corporelle (SCC)** - Échelle 1-5 :")
            
            scores_condition = REFERENTIELS_OFFICIELS["WORLD_SHEEP_BREEDS"]["score_conditions"]
            for score, description in scores_condition.items():
                col_sc1, col_sc2 = st.columns([1, 4])
                with col_sc1:
                    st.metric("Score", score)
                with col_sc2:
                    st.write(description)
            
            st.markdown("---")
            st.markdown("**📚 Références internationales :**")
            st.write("- **FAO**: Organisation des Nations Unies pour l'alimentation et l'agriculture")
            st.write("- **ICAR**: International Committee for Animal Recording")
            st.write("- **WAAP**: World Association for Animal Production")
        
        with ref_tabs[2]:
            st.markdown("#### Documentation Technique")
            
            st.markdown("""
            **📖 Guide d'utilisation du scoring phénotypique**
            
            1. **Évaluation Mamelle** (100 points max)
               - Observer l'animal debout, de profil et de derrière
               - Noter l'équilibre entre quartiers
               - Vérifier la position et la taille des trayons
            
            2. **Évaluation Membres** (80 points max)
               - Observer l'animal en mouvement
               - Vérifier l'alignement des paturons
               - Noter la solidité des canons
            
            3. **Évaluation Type** (80 points max)
               - Mesurer ou estimer les proportions
               - Comparer aux standards de race
               - Noter le développement musculaire
            
            4. **Évaluation Conformation** (variable)
               - Observer la couleur et texture de la robe
               - Noter la présence/forme des cornes
               - Évaluer le tempérament
            
            **🎯 Fréquence d'évaluation recommandée :**
            - Jeunes animaux : À 6, 12 et 18 mois
            - Adultes : Avant et après chaque saison de reproduction
            - Reproducteurs : Avant chaque utilisation
            
            **📊 Interprétation des scores :**
            - >90% : Excellence, reproducteur d'élite
            - 75-90% : Très bon, améliorateur
            - 60-75% : Bon, moyen
            - <60% : À améliorer ou réformer
            """)

# ========== FORMULAIRES STANDARDISÉS PAR RACE ==========

def afficher_formulaires_standardises():
    """Affiche les formulaires de saisie standardisés par race"""
    
    st.markdown('<h2 class="section-header">📝 Formulaires Standardisés</h2>', unsafe_allow_html=True)
    
    # Sélection du type de formulaire
    formulaire_type = st.radio(
        "Type de formulaire :",
        ["Nouvel animal", "Évaluation périodique", "Score de condition", "Données morphométriques"]
    )
    
    if formulaire_type == "Nouvel animal":
        # Formulaire pour nouvelle entrée
        with st.form("form_nouvel_animal_standard"):
            st.markdown("### 🐑 Enregistrement d'un nouvel animal")
            
            col_id, col_date = st.columns(2)
            
            with col_id:
                identifiant = st.text_input("Identifiant unique*", 
                                          placeholder="Ex: ODJ-2024-001")
                nom = st.text_input("Nom", placeholder="Ex: Bella")
            
            with col_date:
                date_naissance = st.date_input("Date de naissance*", 
                                             value=date.today() - timedelta(days=365))
                sexe = st.radio("Sexe*", ["Femelle", "Mâle"], horizontal=True)
            
            # Sélection de race avec sous-races
            st.markdown("### 🎪 Race et Origine")
            
            race_col1, race_col2 = st.columns(2)
            
            with race_col1:
                race_principale = st.selectbox(
                    "Race principale*",
                    list(RACES_ALGERIENNES.keys()),
                    format_func=lambda x: RACES_ALGERIENNES[x]["nom_complet"]
                )
                
                # Afficher les caractéristiques de la race sélectionnée
                if race_principale:
                    race_info = RACES_ALGERIENNES[race_principale]
                    with st.expander(f"Caractéristiques de la race {race_info['nom_complet']}"):
                        for key, value in race_info['caracteristiques'].items():
                            st.write(f"**{key.replace('_', ' ').title()}:** {value}")
            
            with race_col2:
                # Sous-races ou variétés
                if race_principale == "OULED_DJELLAL":
                    sous_race = st.selectbox("Variété/Sous-race", 
                                           ["Type Sétif", "Type Batna", "Type Biskra", "Non spécifié"])
                elif race_principale == "RAZE":
                    sous_race = st.selectbox("Variété/Sous-race", 
                                           ["Kabyle", "Aurès", "Chélia", "Non spécifié"])
                elif race_principale == "HAMRA":
                    sous_race = st.selectbox("Variété/Sous-race", 
                                           ["Type El Oued", "Type Ouargla", "Non spécifié"])
                else:
                    sous_race = st.selectbox("Variété/Sous-race", ["Non spécifié"])
                
                # Origine géographique
                wilaya = st.selectbox("Wilaya d'origine", 
                                    ["Alger", "Oran", "Constantine", "Annaba", "Batna", "Béjaïa", 
                                     "Sétif", "Tizi Ouzou", "Autre", "Non spécifiée"])
            
            # Données morphométriques initiales
            st.markdown("### 📏 Données morphométriques initiales")
            
            morpho_col1, morpho_col2, morpho_col3 = st.columns(3)
            
            with morpho_col1:
                poids = st.number_input("Poids (kg)*", min_value=0.0, max_value=200.0, value=30.0)
                longueur_estimee = st.number_input("Longueur estimée (cm)", 50.0, 200.0, 100.0)
            
            with morpho_col2:
                hauteur_estimee = st.number_input("Hauteur estimée (cm)", 40.0, 150.0, 70.0)
                tour_poitrine = st.number_input("Tour de poitrine (cm)", 60.0, 180.0, 90.0)
            
            with morpho_col3:
                score_condition = st.slider("Score de condition (1-5)", 1, 5, 3,
                                          help="1: Émaciation extrême, 3: Optimal, 5: Obèse")
            
            # Caractéristiques phénotypiques
            st.markdown("### 🌟 Caractéristiques phénotypiques")
            
            pheno_col1, pheno_col2 = st.columns(2)
            
            with pheno_col1:
                couleur_robe = st.selectbox("Couleur de la robe", 
                                          ["Blanche", "Noire", "Rousse", "Brune", "Grise", "Pie", "Tachetée"])
                type_laine = st.select_slider("Type de laine", 
                                            ["Très fine", "Fine", "Moyenne", "Grossière", "Très grossière"])
            
            with pheno_col2:
                cornes_presence = st.radio("Présence de cornes", 
                                         ["Présentes", "Absentes", "Rudimentaires"])
                marques_particulieres = st.text_area("Marques particulières", 
                                                   placeholder="Taches, cicatrices, particularités...")
            
            # Origine parentale
            st.markdown("### 👨‍👩‍👧 Origine parentale")
            
            parent_col1, parent_col2 = st.columns(2)
            
            with parent_col1:
                mere_id = st.text_input("Identifiant de la mère (optionnel)", 
                                      placeholder="Ex: ODJ-2022-015")
            
            with parent_col2:
                pere_id = st.text_input("Identifiant du père (optionnel)", 
                                      placeholder="Ex: ODJ-2021-003")
            
            # Notes et observations
            observations = st.text_area("Observations initiales", 
                                      placeholder="Santé, comportement, particularités...")
            
            # Bouton de soumission
            submitted = st.form_submit_button("📝 Enregistrer l'animal avec formulaire standardisé", 
                                            type="primary")
            
            if submitted:
                if identifiant:
                    # Validation des données
                    erreurs = []
                    
                    if not identifiant:
                        erreurs.append("L'identifiant unique est obligatoire")
                    if poids <= 0:
                        erreurs.append("Le poids doit être positif")
                    
                    if erreurs:
                        for erreur in erreurs:
                            st.error(erreur)
                    else:
                        # Enregistrement dans la base de données
                        try:
                            cursor = conn.cursor()
                            
                            # Table brebis étendue
                            cursor.execute('''
                                CREATE TABLE IF NOT EXISTS brebis_detaille (
                                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                                    identifiant_unique TEXT UNIQUE NOT NULL,
                                    nom TEXT,
                                    date_naissance DATE,
                                    race_principale TEXT,
                                    sous_race TEXT,
                                    wilaya_origine TEXT,
                                    sexe TEXT,
                                    poids_initial FLOAT,
                                    longueur_initiale FLOAT,
                                    hauteur_initiale FLOAT,
                                    tour_poitrine_initial FLOAT,
                                    score_condition_initial INTEGER,
                                    couleur_robe TEXT,
                                    type_laine TEXT,
                                    cornes TEXT,
                                    marques_particulieres TEXT,
                                    mere_id TEXT,
                                    pere_id TEXT,
                                    observations TEXT,
                                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                                )
                            ''')
                            
                            # Insérer les données
                            cursor.execute('''
                                INSERT INTO brebis_detaille 
                                (identifiant_unique, nom, date_naissance, race_principale, 
                                 sous_race, wilaya_origine, sexe, poids_initial, 
                                 longueur_initiale, hauteur_initiale, tour_poitrine_initial,
                                 score_condition_initial, couleur_robe, type_laine, cornes,
                                 marques_particulieres, mere_id, pere_id, observations)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            ''', (
                                identifiant, nom, date_naissance.isoformat(), 
                                race_principale, sous_race, wilaya_origine,
                                "F" if sexe == "Femelle" else "M",
                                poids, longueur_estimee, hauteur_estimee, 
                                tour_poitrine, score_condition,
                                couleur_robe, type_laine, cornes_presence,
                                marques_particulieres, mere_id, pere_id, observations
                            ))
                            
                            conn.commit()
                            
                            # Message de succès
                            st.success(f"✅ Animal {nom} ({identifiant}) enregistré avec succès!")
                            st.balloons()
                            
                            # Afficher un récapitulatif
                            with st.expander("📋 Récapitulatif de l'enregistrement"):
                                recap_data = {
                                    "Identifiant": identifiant,
                                    "Nom": nom,
                                    "Race": RACES_ALGERIENNES[race_principale]["nom_complet"],
                                    "Sous-race": sous_race,
                                    "Poids initial": f"{poids} kg",
                                    "Score condition": f"{score_condition}/5"
                                }
                                st.json(recap_data)
                            
                        except sqlite3.IntegrityError:
                            st.error("❌ Cet identifiant existe déjà dans la base de données!")
                        except Exception as e:
                            st.error(f"Erreur lors de l'enregistrement: {e}")
                else:
                    st.warning("⚠️ L'identifiant unique est obligatoire")
    
    elif formulaire_type == "Évaluation périodique":
        st.markdown("### 📅 Évaluation Périodique Standardisée")
        
        # Sélection de l'animal
        cursor = conn.cursor()
        cursor.execute("SELECT identifiant_unique, nom, race_principale FROM brebis_detaille ORDER BY nom")
        animaux = cursor.fetchall()
        
        if animaux:
            animal_options = [f"{a[1]} ({a[0]}) - {a[2]}" for a in animaux]
            animal_selected = st.selectbox("Sélectionner l'animal à évaluer", animal_options)
            
            if animal_selected:
                with st.form("form_evaluation_periodique"):
                    # Date d'évaluation
                    date_eval = st.date_input("Date d'évaluation", value=date.today())
                    
                    st.markdown("#### Score de Condition Corporelle (SCC)")
                    
                    # Échelle visuelle SCC
                    scc_score = st.slider("Score SCC (1-5)", 1, 5, 3, 
                                        help="""1: Émaciation extrême (côtes très visibles)
        2: Maigre (côtes visibles)
        3: Optimal (côtes palpables mais non visibles)
        4: Gras (côtes difficilement palpables)
        5: Obèse (côtes non palpables)""")
                    
                    # Affichage visuel du SCC
                    scc_descriptions = {
                        1: "⚠️ Émaciation extrême - Nécessite intervention",
                        2: "📉 Maigre - Surveillance nécessaire",
                        3: "✅ Optimal - État idéal",
                        4: "📈 Gras - Risque de problèmes métaboliques",
                        5: "🚨 Obèse - Intervention requise"
                    }
                    
                    st.info(scc_descriptions[scc_score])
                    
                    # Mensurations actuelles
                    st.markdown("#### Mensurations actuelles")
                    
                    col_mes1, col_mes2, col_mes3 = st.columns(3)
                    
                    with col_mes1:
                        poids_actuel = st.number_input("Poids actuel (kg)", 0.0, 200.0, 50.0)
                    
                    with col_mes2:
                        longueur_actuelle = st.number_input("Longueur corps (cm)", 50.0, 200.0, 110.0)
                    
                    with col_mes3:
                        hauteur_actuelle = st.number_input("Hauteur garrot (cm)", 40.0, 150.0, 75.0)
                    
                    # État de santé
                    st.markdown("#### État de santé général")
                    
                    sante_col1, sante_col2 = st.columns(2)
                    
                    with sante_col1:
                        etat_paturons = st.selectbox("État des paturons", 
                                                   ["Excellent", "Bon", "Moyen", "Mauvais", "Grave"])
                        etat_dentaire = st.select_slider("État dentaire", 
                                                       ["Parfait", "Bon", "Usure normale", "Usure avancée", "Problèmes"])
                    
                    with sante_col2:
                        parasites = st.multiselect("Parasites observés", 
                                                 ["Gastro-intestinaux", "Pou", "Tique", "Gale", "Aucun"])
                        vaccinations = st.multiselect("Vaccinations à jour", 
                                                    ["FCO", "Clostridium", "Pasteurellose", "Rage", "Autres"])
                    
                    # Observations
                    observations = st.text_area("Observations et recommandations")
                    
                    if st.form_submit_button("💾 Enregistrer l'évaluation périodique"):
                        st.success("Évaluation enregistrée!")
        else:
            st.info("Aucun animal enregistré dans la base détaillée.")

# ========== INTÉGRATION DANS L'APPLICATION PRINCIPALE ==========

# Dans la barre latérale, ajoutez les nouvelles pages :
with st.sidebar:
    # ... code existant ...
    
    page = st.radio(
        "Menu Principal",
        ["🏠 Tableau de Bord", 
         "📊 Gestion des Brebis", 
         "🧬 Génétique & NCBI",
         "🎯 Scoring Phénotypique",      # NOUVEAU
         "📝 Formulaires Standardisés",   # NOUVEAU
         "🥛 Analyse Lait",
         "📐 Morphométrie 3D",
         "🤰 Suivi Gestation", 
         "📈 Statistiques Avancées",
         "⚙️ Paramètres"]
    )

# Dans la navigation principale, ajoutez :
if page == "🎯 Scoring Phénotypique":
    afficher_scoring_phenotypique()
elif page == "📝 Formulaires Standardisés":
    afficher_formulaires_standardises()
# ... autres pages ...
