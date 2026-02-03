"""
Module: Gestion de la base de données
"""

import sqlite3
from datetime import datetime, date
from typing import List, Dict, Optional
import pandas as pd

class DatabaseManager:
    """Gestionnaire de base de données SQLite"""
    
    def __init__(self, db_path: str = "data/ovin_manager.db"):
        self.db_path = db_path
        self.conn = None
        self.connect()
    
    def connect(self):
        """Établit la connexion à la base de données"""
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row
            return True
        except sqlite3.Error as e:
            print(f"Erreur connexion DB: {e}")
            return False
    
    def execute_query(self, query: str, params: tuple = ()):
        """Exécute une requête SQL"""
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            self.conn.commit()
            return cursor
        except sqlite3.Error as e:
            print(f"Erreur requête: {e}")
            return None
    
    def fetch_all(self, query: str, params: tuple = ()):
        """Récupère tous les résultats d'une requête"""
        cursor = self.execute_query(query, params)
        if cursor:
            return cursor.fetchall()
        return []
    
    def fetch_one(self, query: str, params: tuple = ()):
        """Récupère un seul résultat"""
        cursor = self.execute_query(query, params)
        if cursor:
            return cursor.fetchone()
        return None
    
    def close(self):
        """Ferme la connexion"""
        if self.conn:
            self.conn.close()

def init_database(db_manager: DatabaseManager):
    """Initialise la base de données avec les tables"""
    
    tables = [
        # Table brebis
        """
        CREATE TABLE IF NOT EXISTS brebis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            identifiant_unique TEXT UNIQUE NOT NULL,
            nom TEXT,
            date_naissance DATE,
            race TEXT,
            sexe TEXT,
            statut TEXT DEFAULT 'active',
            notes TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        """,
        
        # Table gestations
        """
        CREATE TABLE IF NOT EXISTS gestations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_eponge DATE,
            date_retrait_eponge DATE,
            date_insemination DATE,
            date_mise_bas_prevu DATE,
            date_mise_bas_reel DATE,
            nombre_agneaux_prevus INTEGER DEFAULT 1,
            nombre_agneaux_nes INTEGER DEFAULT 0,
            statut TEXT DEFAULT 'en_cours',
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
        """,
        
        # Table suivi_medical
        """
        CREATE TABLE IF NOT EXISTS suivi_medical (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_intervention DATE,
            type_intervention TEXT,
            produit TEXT,
            dose TEXT,
            veterinaire TEXT,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
        """,
        
        # Table caracteres_morpho
        """
        CREATE TABLE IF NOT EXISTS caracteres_morpho (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            date_mesure DATE,
            caractere TEXT,
            valeur REAL,
            unite TEXT,
            methode_mesure TEXT DEFAULT 'manuel',
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
        """,
        
        # Table sequences_genetiques
        """
        CREATE TABLE IF NOT EXISTS sequences_genetiques (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            brebis_id INTEGER,
            accession_number TEXT UNIQUE,
            sequence_type TEXT,
            longueur INTEGER,
            date_sequencage DATE,
            laboratoire TEXT,
            notes TEXT,
            FOREIGN KEY (brebis_id) REFERENCES brebis (id)
        )
        """
    ]
    
    for table_sql in tables:
        db_manager.execute_query(table_sql)
    
    print("✅ Base de données initialisée")
