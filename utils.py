"""
Module: Utilitaires divers
"""

import pandas as pd
import io
from datetime import datetime

def generer_rapport_pdf(data, title: str) -> io.BytesIO:
    """Génère un rapport PDF (simulé)"""
    # En production, utiliser reportlab ou autre lib
    rapport = f"Rapport: {title}\nDate: {datetime.now()}\n\n"
    
    if isinstance(data, pd.DataFrame):
        rapport += data.to_string()
    else:
        rapport += str(data)
    
    buffer = io.BytesIO(rapport.encode())
    buffer.seek(0)
    return buffer

def exporter_excel(data, nom_feuille: str = "Données") -> io.BytesIO:
    """Exporte des données en Excel"""
    if isinstance(data, pd.DataFrame):
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            data.to_excel(writer, sheet_name=nom_feuille, index=False)
        output.seek(0)
        return output
    else:
        buffer = io.BytesIO()
        buffer.write(b"Donnees Excel")
        buffer.seek(0)
        return buffer
