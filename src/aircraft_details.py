import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Gestion automatique du chemin vers le fichier Excel
base_dir = os.path.dirname(os.path.abspath(__file__))
excel_path = os.path.join(base_dir, "data", "Aircraft_Data.xlsx")

# Lecture
try:
    df = pd.read_excel(excel_path)
except FileNotFoundError:
    print(f"Erreur : Le fichier {excel_path} est introuvable.")
    exit()

var_x = "MTOW"
var_y = "ZFW"

# Correction : On ne garde que les colonnes numériques pour le pairplot
df_numerique = df.select_dtypes(include=['float64', 'int64'])

# Création du Pairplot
sns.pairplot(df_numerique.iloc[:, 0:9]) # On limite aux 9 premières colonnes numériques

# Création du Scatter plot spécifique
plt.figure(figsize=(8, 6))
plt.scatter(df[var_x], df[var_y], c='blue', alpha=0.6)
plt.title(f"Comparaison : {var_x} vs {var_y}")
plt.xlabel(var_x)
plt.ylabel(var_y)
plt.grid(True)

plt.show()