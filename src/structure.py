import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
import os
import pandas as pd # Nécessaire pour lire les données réelles

# --- GESTION DES CHEMINS ---
# Permet de trouver le dossier 'data' même si on lance le script d'ailleurs
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, 'data'))

## Données pour le flambage local
try:
    from Data_k_bh import bh, k
    # On garde les données telles quelles
    bh_data = bh
    k_data = k
except ImportError:
    print("! ATTENTION : Module data.Data_k_bh introuvable. Utilisation de valeurs par défaut.")
    bh_data = [0, 100]
    k_data = [4, 4]

# Création de la fonction d'interpolation
get_k_factor = interp1d(bh_data, k_data, kind='linear', fill_value="extrapolate")

## Données et contraintes
E_ALU = 70e9       
NU_ALU = 0.33      

h_sect = 0.07      # 7 cm
b_sect = 0.03      # 3 cm
e_sect = 0.004     # 4 mm

L_poutre = 20.0    
x_F_pos = 4.0      
F_load = 200.0     
q_load = -15.0     

## ---------------- FONCTIONS DE CALCUL ---------------- ##

def calcul_section(h, b, e):
    """
    Calcule S (Aire), Ix (Inertie flexion), et Iy pour une section en I.
    """
    # Aire : 2*semelles + ame
    S = 2*b*e + h*e
    # Inertie Ix (autour axe neutre horizontal)
    I_x = e*( h**3/12 + b*e**2/6 + b*(h + e)**2/2)
    # Inertie Iy
    I_y = e/6*( h*e**2/2 + b**3)
    return S, I_x, I_y

def calcul_poutre_console_superposition(x_array, L, a, F, q, E, I):
    """
    Calcule M(x) et y(x) par superposition.
    """
    y_tot = np.zeros_like(x_array)
    M_tot = np.zeros_like(x_array)

    # --- Cas 1 : Charge répartie q ---
    M_q = -(q / 2) * (L - x_array)**2
    # Formule de la flèche pour charge répartie (attention au facteur L parfois oublié)
    y_q = (q / (24 * E * I)) * x_array**2 * (x_array**2 - 4*L*x_array + 6*L**2)

    # --- Cas 2 : Force ponctuelle F ---
    M_F = np.zeros_like(x_array)
    y_F = np.zeros_like(x_array)

    for i, x in enumerate(x_array):
        if x <= a:
            M_F[i] = -F * (a - x)
            y_F[i] = (F / (6 * E * I)) * x**2 * (3*a - x)
        else:
            M_F[i] = 0
            y_F[i] = (F / (6 * E * I)) * a**2 * (3*x - a)

    return M_q + M_F, y_q + y_F

def calcul_force_critique(L, E, I, S, nu, b, h, e):
    """
    Calcule la force critique (Euler vs Johnson/Crippling).
    """
    # 1. Flambage Global (Euler)
    K = 0.5
    L_eq = K * L
    if L_eq <= 0: L_eq = 0.001
    F_euler = (np.pi**2 * E * I) / (L_eq**2)

    # 2. Voilement Local (Crippling)
    ratio_bh = b / h
    k_crippling = float(get_k_factor(ratio_bh))
    
    sigma_cc = ((np.pi * e) / h)**2 * (E / (12 * (1 - nu**2))) * k_crippling

    # 3. Flambage Réel (Johnson)
    rho = np.sqrt(I / S)
    elancement = L_eq / rho
    
    # Sécurité mathématique
    if sigma_cc <= 0:
        elancement_critique = 1e9
    else:
        elancement_critique = np.pi * np.sqrt(2 * E / sigma_cc)

    if elancement >= elancement_critique:
        return F_euler, "Euler (Global)", elancement
    else:
        sigma_johnson = sigma_cc - (sigma_cc**2 * elancement**2) / (4 * np.pi**2 * E)
        return sigma_johnson * S, "Johnson (Local/Plastique)", elancement

def optimisation_section(L, F_cible, E, nu):
    """
    Cherche la section la plus légère qui résiste à la charge F_cible.
    """
    solutions = []
    # e entre 1 mm et 10 mm
    for e in range(1, 10, 1):
        e_eval = e/(10**3)
        # b entre 1 cm et 20 cm
        for b in range (1, 20, 1):
            b_eval = b/(10**2)
            # h entre 1 cm et 20 cm
            for h in range(1, 20, 1):
                h_eval = h/(10**2)
                
                S, Ix, Iy = calcul_section(h_eval, b_eval, e_eval)
                f_max = calcul_force_critique(L, E, Ix, S, nu, b_eval, h_eval, e_eval)[0]
                
                # CONDITION CRITIQUE : La poutre doit résister (F_max >= Charge)
                if f_max >= F_cible:
                    solutions.append({'S': S, 'e': e_eval, 'b': b_eval, 'h': h_eval, 'F_max': f_max})

    if not solutions:
        return "Aucune section trouvée."
    
    # On trie par surface croissante (le plus léger en premier)
    best = sorted(solutions, key=lambda x: x['S'])[0]
    return (f"OPTIMUM : e={best['e']*1000}mm, b={best['b']*100}cm, h={best['h']*100}cm "
            f"(Masse lin.={best['S']*2700:.2f} kg/m, F_crit={best['F_max']:.1f}N)")

def calcul_masse_aile_leclerc(MTOW, S, b, epsilon, phi_25_rad, e_r):
    """
    Calcule la masse de l’aile selon la formule de F. Leclerc (2002).
    """
    # Allongement
    if S <= 0: return 0
    lamb = b**2 / S

    # Sécurité cosinus
    cos_phi = np.cos(phi_25_rad)
    if abs(cos_phi) < 0.001: cos_phi = 0.001
        
    term_1 = (MTOW * lamb * (epsilon + 1)) / (e_r * cos_phi)

    # Calcul final (Masse Wing)
    Mw = 0.197 * (term_1**0.6) * (S**0.3)

    return Mw


## ---------------- EXÉCUTION ET AFFICHAGE ---------------- ##

print("RAPPORT DE LA STRUCTURE AVION \n")

# PARTIE 1 : SECTION
S, Ix, Iy = calcul_section(h_sect, b_sect, e_sect)
print(f"1. PROPRIÉTÉS DE LA SECTION")
print(f"   - Aire (S) : {S:.6e} m²")
print(f"   - Ix : {Ix:.6e} m^4")
print(f"   - Iy : {Iy:.6e} m^4")
print("-" * 50)

# PARTIE 2 : FLEXION
x_vals = np.linspace(0, L_poutre, 200)
M_res, y_res = calcul_poutre_console_superposition(x_vals, L_poutre, x_F_pos, F_load, q_load, E_ALU, Ix)

idx_max_M = np.argmax(np.abs(M_res))
idx_max_y = np.argmax(np.abs(y_res))

print(f"2. RÉSULTATS FLEXION")
print(f"   - Moment Max : {abs(M_res[idx_max_M]):.2f} N.m")
print(f"   - Flèche Max : {abs(y_res[idx_max_y])*1000:.2f} mm")
print("-" * 50)

# PARTIE 3 : FLAMBAGE
print(f"3. ANALYSE DU FLAMBAGE")
for L_test in [0.5, 2.0]:
    F_crit, mode_ruine, lam = calcul_force_critique(L_test, E_ALU, Ix, S, NU_ALU, b_sect, h_sect, e_sect)
    print(f"   L={L_test}m : F_crit={F_crit:.1f} N ({mode_ruine})")

# PARTIE 4 : OPTIMISATION
print("-" * 50)
print("4. RÉSULTAT OPTIMISATION (Ex 5)")
print(optimisation_section(4, 100, E_ALU, NU_ALU))

# PARTIE 5 : CALCUL MASSE AILE (DONNÉES RÉELLES)
print("-" * 50)
print("5. ESTIMATION MASSE AILE (Données: Aircraft_Data.xlsx)")

try:
    # Lecture du fichier Excel
    file_path = os.path.join(current_dir, "data", "Aircraft_Data.xlsx")
    
    # On gère le cas où le fichier serait un CSV (fallback) ou un Excel
    if os.path.exists(file_path):
        df_avions = pd.read_excel(file_path)
    else:
        # Essai avec extension csv si xlsx manquant (cas fréquents lors des tests)
        df_avions = pd.read_csv(os.path.join(current_dir, "data", "Aircraft_Data.csv"))

    # --- SÉLECTION DE L'AVION (Ligne 0 : A300 B4) ---
    avion = df_avions.iloc[0] 
    
    # --- EXTRACTION DES PARAMÈTRES ---
    # Note : Les colonnes Surface et Envergure n'ont pas de nom dans l'entête du fichier fourni,
    # on utilise donc leur position (index 13 et 14 vérifiés).
    
    mtow_val = avion['MTOW']       # Masse Max (kg)
    surface_val = avion.iloc[13]   # Surface S (m²) - Colonne index 13
    envergure_val = avion.iloc[14] # Envergure b (m) - Colonne index 14
    
    # Angle de flèche (On prend la colonne 'Avant' en degrés)
    phi_deg = avion['Avant']  
    phi_rad = np.radians(phi_deg)
    
    # Estimation de l'épaisseur à l'emplanture (e_r)
    # On approxime e_r via la Corde à l'emplanture ('Emp') * ratio t/c (ex: 15%)
    corde_emplanture = avion['Emp'] 
    e_r_val = 0.15 * corde_emplanture 
    
    epsilon_val = 0.0 # Aile sèche par défaut

    # --- CALCUL ---
    masse_estimee = calcul_masse_aile_leclerc(mtow_val, surface_val, envergure_val, epsilon_val, phi_rad, e_r_val)
    
    # --- AFFICHAGE STRICT (UNIQUEMENT LA MASSE) ---
    print(f"Masse de l'aile calculée : {masse_estimee:.2f} kg")

except Exception as e:
    print(f"Impossible de lire les données avion : {e}")

# --- GRAPHIQUES ---
l_range = np.linspace(0.1, 5.0, 100)
forces_critiques = []
elancements = []

ratio = b_sect / h_sect
k_val = float(get_k_factor(ratio))
sigma_cc_temp = ((np.pi * e_sect) / h_sect)**2 * (E_ALU / (12 * (1 - NU_ALU**2))) * k_val
lambda_crit = np.pi * np.sqrt(2 * E_ALU / sigma_cc_temp) if sigma_cc_temp > 0 else 0

for l in l_range:
    f_cr, _, lam = calcul_force_critique(l, E_ALU, Ix, S, NU_ALU, b_sect, h_sect, e_sect)
    forces_critiques.append(f_cr)
    elancements.append(lam)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), constrained_layout=True)

ax1.plot(x_vals, M_res, 'r-')
ax1.set_title('Moment Fléchissant (N.m)')
ax1.grid(True)

ax2.plot(x_vals, y_res * 1000, 'b-')
ax2.set_title('Déformée (mm)')
ax2.invert_yaxis()
ax2.grid(True)

ax3.plot(elancements, forces_critiques, 'g-', label='Force Critique')
if lambda_crit > 0:
    ax3.axvline(x=lambda_crit, color='k', linestyle='--', label='Transition')
ax3.set_title("Flambage : Force Critique")
ax3.set_xlabel("Élancement")
ax3.legend()
ax3.grid(True)

plt.show()