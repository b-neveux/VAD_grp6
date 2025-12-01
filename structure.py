import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

##Données et contraintes

# Matériau : Alliage Aluminium (Valeurs standards)
E_ALU = 70e9       # 70 GPa
NU_ALU = 0.33      # Coefficient de Poisson

# Géométrie de la section (Exercice 3 & 4 : Section en I)
h_sect = 0.07      # 7 cm = 0.07 m
b_sect = 0.03      # 3 cm = 0.03 m
e_sect = 0.004     # 4 mm = 0.004 m

# Chargement (Exercice 3 - Cas c/)
L_poutre = 20.0    # Longueur totale (m)
x_F_pos = 4.0      # Position de la force ponctuelle (a=4m)
F_load = 200.0     # Force ponctuelle (N) - vers le bas
q_load = -15.0     # Force linéique (N/m) - vers le haut (donc négatif)

##Données pour le flambage local
# Ces listes permettent de calculer le coefficient k pour le "crippling" (voilement)
from Sources.Data_k_bh import bh, k

bh_data = bh
k_data = k

# Création de la fonction d'interpolation linéaire
get_k_factor = interp1d(bh_data, k_data, kind='linear', fill_value="extrapolate")

##Fonctions de calcul

def calcul_section(h, b, e):
    """
    Calcule S (Aire), Ix (Inertie flexion), et Rho (Rayon giration) pour une section en I.
    """
    # Aire : 2*petit rectangle + grand rectangle
    S = 2*b*e + h*e

    # Inertie Ix (autour axe neutre horizontal)
    I_x = e*( h**3/12 + b*e**2/6 + b*(h + e)**2/2)

    I_y=e/6*( h*e**2/2 + b**3)

    return S, I_x, I_y

def calcul_poutre_console_superposition(x_array, L, a, F, q, E, I):
    """
    Calcule M(x) et y(x) par superposition des cas C01 (répartie) et C02 (ponctuelle).
    """
    # Initialisation des vecteurs
    y_tot = np.zeros_like(x_array)
    M_tot = np.zeros_like(x_array)

    # --- Cas 1 : Charge répartie q (C01) ---
    # Formules du cours (avec Q = q*L)
    # M(x) = -q/2 * (L-x)^2
    # y(x) = (q / 24EI) * x^2 * (x^2 - 4Lx + 6L^2)
    M_q = -(q / 2) * (L - x_array)**2
    y_q = (q / (24 * E * I *L)) * x_array**2 * (x_array**2 - 4*L*x_array + 6*L**2)

    # --- Cas 2 : Force ponctuelle F en a (C02) ---
    M_F = np.zeros_like(x_array)
    y_F = np.zeros_like(x_array)

    # Formules par morceaux
    for i, x in enumerate(x_array):
        if x <= a:
            M_F[i] = -F * (a - x)
            y_F[i] = (F / (6 * E * I)) * x**2 * (3*a - x)
        else:
            M_F[i] = 0
            y_F[i] = (F / (6 * E * I)) * a**2 * (3*x - a)

    # Superposition
    y_tot = y_q + y_F
    M_tot = M_q + M_F

    return M_tot, y_tot

def calcul_force_critique(L, E, I, S, nu, b, h, e):
    """
    Calcule la force critique (Euler vs Johnson/Crippling).
    """
    # 1. Flambage Global (Euler) - Hypothèse Bi-encastré K=0.5
    K = 0.5
    L_eq = K * L
    F_euler = (np.pi**2 * E * I) / (L_eq**2)

    # 2. Voilement Local (Crippling)
    # On récupère k via interpolation des données fournies
    ratio_bh = b / h
    k_crippling = float(get_k_factor(ratio_bh))

    # Contrainte critique de voilement (Formule cours)
    # Sigma_cc = (pi * e / h)^2 * E / (12 * (1 - nu^2)) * k
    sigma_cc = ((np.pi * e) / h)**2 * (E / (12 * (1 - nu**2))) * k_crippling

    # 3. Flambage Réel (Interaction Johnson-Euler)
    rho = np.sqrt(I / S)
    elancement = L_eq / rho
    # Limite d'élancement entre plastique et élastique
    elancement_critique = np.pi * np.sqrt(2 * E / sigma_cc)

    F_critique_finale = 0
    mode = ""

    if elancement >= elancement_critique:
        # Domaine élastique -> Euler
        F_critique_finale = F_euler
        mode = "Euler (Global)"
    else:
        # Domaine plastique/local -> Johnson
        # Sigma_johnson = Sigma_cc - (Sigma_cc^2 * lambda^2) / (4 * pi^2 * E)
        sigma_johnson = sigma_cc - (sigma_cc**2 * elancement**2) / (4 * np.pi**2 * E)
        F_critique_finale = sigma_johnson * S
        mode = "Johnson (Local/Plastique)"

    return F_critique_finale, mode, elancement

##Exécution et affichage

print("RAPPORT DE LA STRUCTURE AVION \n")

#PARTIE 1 : SECTION (EXERCICE 1)
S, Ix, Iy = calcul_section(h_sect, b_sect, e_sect)
print(f"1. PROPRIÉTÉS DE LA SECTION (I de {h_sect*100}x{b_sect*100} cm, ép {e_sect*1000} mm)")
print(f"   - Aire (S) : {S:.6e} m²")
print(f"   - Moment Quadratique (Ix) : {Ix:.6e} m^4")
print(f"   - Moment Quadratique (Iy) : {Iy:.6e} m^4")
print("-" * 50)

#PARTIE 2 : FLEXION (EXERCICES 2 & 3)
# Discrétisation
x_vals = np.linspace(0, L_poutre, 200)
M_res, y_res = calcul_poutre_console_superposition(x_vals, L_poutre, x_F_pos, F_load, q_load, E_ALU, Ix)

# Résultats max
idx_max_M = np.argmax(np.abs(M_res))
idx_max_y = np.argmax(np.abs(y_res))

print(f"2. RÉSULTATS FLEXION (L={L_poutre}m, F={F_load}N à {x_F_pos}m, q={q_load}N/m)")
print(f"   - Moment Fléchissant Max : {abs(M_res[idx_max_M]):.2f} N.m (à x={x_vals[idx_max_M]:.2f} m)")
print(f"   - Flèche Maximale (bout) : {abs(y_res[idx_max_y])*1000:.2f} mm (à x={x_vals[idx_max_y]:.2f} m)")
print("-" * 50)

# Graphique Flexion
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
# Moment
ax1.plot(x_vals, M_res, 'r-', linewidth=2)
ax1.fill_between(x_vals, M_res, 0, color='r', alpha=0.1)
ax1.set_title('Diagramme du Moment Fléchissant $M_f(x)$')
ax1.set_ylabel('Moment (N.m)')
ax1.grid(True)
# Déformée
ax2.plot(x_vals, y_res * 1000, 'b-', linewidth=2) # Conversion en mm
ax2.set_title('Déformée de la poutre $y(x)$')
ax2.set_xlabel('Position x (m)')
ax2.set_ylabel('Flèche (mm)')
ax2.invert_yaxis() # Flèche positive vers le bas
ax2.grid(True)
plt.tight_layout()
plt.show()

#PARTIE 3 : FLAMBAGE (EXERCICE 4)
print(f"3. ANALYSE DU FLAMBAGE (Euler vs Crippling vs Johnson)")
longueurs_test = [0.5, 2.0]

for L_test in longueurs_test:
    F_crit, mode_ruine, lam = calcul_force_critique(L_test, E_ALU, Ix, S, NU_ALU, b_sect, h_sect, e_sect)
    print(f"   -> Pour L = {L_test} m :")
    print(f"      Force Critique = {F_crit:.2f} N")
    print(f"      Élancement = {lam:.2f}")
    print(f"      Mode de ruine dominant : {mode_ruine}")

#PARTIE 4 : GRAPHIQUE ÉVOLUTION FLAMBAGE
# Calcul sur une plage de longueurs
l_range = np.linspace(0.1, 5.0, 100)
forces_critiques = []
elancements = []

# Calcul de la transition théorique pour l'affichage vertical
ratio = b_sect / h_sect
k_val = float(get_k_factor(ratio))
sigma_cc_temp = ((np.pi * e_sect) / h_sect)**2 * (E_ALU / (12 * (1 - NU_ALU**2))) * k_val
lambda_crit = np.pi * np.sqrt(2 * E_ALU / sigma_cc_temp)

for l in l_range:
    f_cr, _, lam = calcul_force_critique(l, E_ALU, Ix, S, NU_ALU, b_sect, h_sect, e_sect)
    forces_critiques.append(f_cr)
    elancements.append(lam)

plt.figure(figsize=(10, 6))
plt.plot(elancements, forces_critiques, 'g-', linewidth=2, label='Force Critique Réelle')
plt.axvline(x=lambda_crit, color='k', linestyle='--', label=rf'Transition Johnson/Euler ($\lambda$={lambda_crit:.1f})')

plt.title("Evolution de la Force Critique de Flambage en fonction de l'élancement")
plt.xlabel(r"Élancement $\lambda$")
plt.ylabel("Force Critique (N)")
plt.legend()
plt.grid(True)
plt.show()