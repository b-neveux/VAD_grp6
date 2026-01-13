import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns
import os

def get_data_path():
    """Récupère le chemin du fichier Excel de manière robuste"""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_dir, "data", "Aircraft_Data.xlsx")

def effectuer_regression_puissance(afficher_graphe=True):
    """
    Régression M_wing = a * MTOW ^ b (inchangé)
    """
    if afficher_graphe:
        print("--- CALCUL DE L'ÉQUATION DE CRÉATION (LOI DE PUISSANCE) ---\n")

    excel_path = get_data_path()
    
    if not os.path.exists(excel_path):
        print(f"ERREUR : Fichier introuvable : {excel_path}")
        return 0, 0

    try:
        df = pd.read_excel(excel_path)
        data_clean = df[['MTOW', 'Wing']].dropna()
        data_clean = data_clean[(data_clean['MTOW'] > 0) & (data_clean['Wing'] > 0)]
        
        X_1d = data_clean['MTOW'].values        
        y_1d = data_clean['Wing'].values        
        
        X_log = np.log(data_clean[['MTOW']].values)
        y_log = np.log(y_1d)
        
        model = LinearRegression()
        model.fit(X_log, y_log)

        b_exponent = model.coef_[0]       
        ln_a = model.intercept_           
        a_coeff = np.exp(ln_a)            

        if afficher_graphe:
            score_r2 = model.score(X_log, y_log)
            print(f"Equation : M_wing = {a_coeff:.4f} * MTOW ^ {b_exponent:.4f}")
            print(f"R² : {score_r2:.4f}")
            
            plt.figure(figsize=(8, 5))
            plt.scatter(X_1d, y_1d, alpha=0.5, label='Données')
            X_plot = np.linspace(X_1d.min(), X_1d.max(), 100)
            plt.plot(X_plot, a_coeff * (X_plot ** b_exponent), 'r', label='Modèle')
            plt.title("Régression Puissance")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()

        return a_coeff, b_exponent

    except Exception as e:
        print(f"Erreur Régression : {e}")
        return 0, 0

def afficher_correlations_mtow():
    """
    Affiche un graphique en barres classant les variables selon leur lien avec la MTOW.
    """
    print("\n--- ANALYSE D'INFLUENCE SUR LA MTOW ---")
    
    excel_path = get_data_path()
    if not os.path.exists(excel_path):
        return

    try:
        df = pd.read_excel(excel_path)
        
        # 1. Garder uniquement le numérique
        df_numeric = df.select_dtypes(include=[np.number]).dropna()
        
        if 'MTOW' not in df_numeric.columns:
            print("Erreur : La colonne 'MTOW' n'est pas présente dans les données numériques.")
            return

        # 2. Calculer la corrélation de toutes les colonnes par rapport à la MTOW
        # On utilise corrwith ou simplement corr()['MTOW']
        correlations = df_numeric.corr()['MTOW']
        
        # On retire la MTOW elle-même (car corrélation de 1 avec elle-même, inutile)
        correlations = correlations.drop('MTOW')
        
        # 3. Trier les résultats par ordre décroissant (du plus corrélé au moins corrélé)
        correlations = correlations.sort_values(ascending=False)
        
        print("Corrélations calculées :")
        print(correlations)

        # 4. Affichage Graphique (Bar Chart Horizontal)
        plt.figure(figsize=(10, 8))
        
        # Création du barplot
        # On utilise une palette de couleurs qui change selon la valeur (bleu négatif, rouge positif)
        sns.barplot(x=correlations.values, y=correlations.index, palette='coolwarm')
        
        plt.title("Corrélation des variables avec la MTOW\n(Plus la barre est longue, plus le lien est fort)")
        plt.xlabel("Coefficient de corrélation (de -1 à 1)")
        plt.ylabel("Variables")
        
        # Ligne verticale à 0 pour bien séparer positif/négatif
        plt.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
        
        plt.grid(axis='x', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Erreur lors de l'analyse MTOW : {e}")

if __name__ == "__main__":
    # 1. Régression classique
    effectuer_regression_puissance(afficher_graphe=True)
    
    # 2. Analyse spécifique : Qui prédit la MTOW ?
    afficher_correlations_mtow()