import numpy as np
import pandas as pd
import sys
import os
from scipy.optimize import minimize, fsolve
from sklearn.linear_model import LinearRegression

# 1. PARAMÈTRES MATÉRIAUX ET CONSTANTES
E_MODULUS = 71.7e9       # [Pa] Alu 7075/2024
YIELD_STRESS = 450e6     # [Pa] Limite élastique
DENSITY = 2800.0         # [kg/m^3] Masse volumique
SAFETY_FACTOR = 1.5      # Facteur de sécurité
POISSON_RATIO = 0.33   

# 2. PARTIE A : FORMULES DE MASSE (THÈSE & RÉGRESSION)
def get_data_path():
    """Récupère le chemin du fichier Excel de manière robuste"""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_dir, "data", "Aircraft_Data.xlsx")

def calcul_regression_custom(mtow_target):
    excel_path = get_data_path()
    if not os.path.exists(excel_path):
        return None

    try:
        df = pd.read_excel(excel_path)
        data_clean = df[['MTOW', 'Wing']].dropna()
        data_clean = data_clean[(data_clean['MTOW'] > 0) & (data_clean['Wing'] > 0)]
        
        X_log = np.log(data_clean[['MTOW']].values)
        y_log = np.log(data_clean['Wing'].values)
        
        model = LinearRegression()
        model.fit(X_log, y_log)
        
        b_exp = model.coef_[0]
        a_coeff = np.exp(model.intercept_)
        
        mass_est = a_coeff * (mtow_target ** b_exp)
        return mass_est, a_coeff, b_exp
    except Exception:
        return None

def calc_phi_50(phi_25_deg, AR, taper):
    """ Convertit la flèche à 25% en flèche à 50% """
    phi_25 = np.radians(phi_25_deg)
    tan_phi_50 = np.tan(phi_25) - (4/AR) * 0.25 * ((1 - taper) / (1 + taper))
    return np.degrees(np.arctan(tan_phi_50))

def formule_these_kroo(MTOW, MZFW, S, b, tc, taper, phi_25, nz=3.75):
    """ Formule de I. Kroo (2001) citée dans la Thèse Roux (p.84). """
    phi_50_deg = calc_phi_50(phi_25, (b**2)/S, taper)
    phi_rad = np.radians(phi_50_deg)
    
    term1 = 20.6 * S
    term2_num = 5.387e-6 * nz * (b**3) * np.sqrt(MTOW * MZFW) * (1 + 2*taper)
    term2_den = tc * (np.cos(phi_rad)**2) * S * (1 + taper)
    
    return term1 + (term2_num / term2_den)

def formule_these_torenbeek(MTOW, MZFW, S, b, tc, taper, phi_25, nz=3.75):
    """ Formule d'E. Torenbeek (1986) citée dans la Thèse Roux (p.86). """
    # Conversion en Impérial pour la formule empirique
    mtow_lbs = MTOW * 2.20462
    mzfw_lbs = MZFW * 2.20462
    b_ft = b * 3.28084
    
    phi_50_deg = calc_phi_50(phi_25, (b**2)/S, taper)
    phi_50 = np.radians(phi_50_deg)
    
    k_no = 1.0 + np.sqrt(1.905 / b_ft)
    k_lambda = (1 + taper)**0.4
    
    # Constante A pour M_WSP = A * (MZFW - Mw)^0.55
    A = 4.58e-3 * k_no * k_lambda * (nz)**0.55 * (b_ft**1.675) * (tc**0.45) * (np.cos(phi_50)**1.325)
    
    def func_to_solve(mw_total_guess_lbs):
        relief_lbs = max(0, mzfw_lbs - mw_total_guess_lbs)
        mw_sp_calc_lbs = A * (relief_lbs ** 0.55)
        # Ratio Structure Primaire / Totale (~0.72)
        mw_total_calc = mw_sp_calc_lbs / 0.72 
        return mw_total_guess_lbs - mw_total_calc

    try:
        mw_final_lbs = fsolve(func_to_solve, 0.1 * mtow_lbs)[0]
    except:
        mw_final_lbs = 0.0
    
    return mw_final_lbs * 0.453592

# 3. PARTIE B : MÉCANIQUE & OPTIMISATION (RDM)
def calculate_bending_moment_root(MTOW, b, nz):
    Lift_total = MTOW * 9.81 * nz
    L_half = Lift_total / 2
    y_cp = (2 * b) / (3 * np.pi) 
    return L_half * y_cp

def check_structure_constraints_scaled(x_scaled, params):
    """
    x_scaled = [t_mm, h_cm, w_cm, n_str_float]
    IMPORTANT: Tout reste en float pour les gradients.
    """
    t_mm, h_cm, w_cm, n_str = x_scaled
    
    # Conversion unités physiques (SI)
    t_skin = t_mm / 1000.0
    h_str  = h_cm / 100.0
    w_str  = w_cm / 100.0
    # On garde n_str en float pour l'optimiseur (approximation continue du nb de lisses)
    
    M_f = params['Moment_Flexion']
    chord = params['chord_root']
    
    h_box = chord * params['tc_ratio'] * 0.9 
    w_box = chord * 0.5                      
    
    # 1. Inertie
    Area_skin = w_box * t_skin
    Area_str_tot = n_str * (h_str * w_str)
    Area_tot = Area_skin + Area_str_tot
    
    # Ixx approx (Caisson rectangulaire + lisses)
    Ixx = 2 * (Area_tot * (h_box / 2)**2)
    
    # 2. Contrainte Flexion
    sigma_f = (M_f * (h_box / 2)) / Ixx
    
    # 3. Marges
    # a) Yield
    margin_yield = (YIELD_STRESS / (SAFETY_FACTOR * sigma_f)) - 1.0
    
    # b) Buckling Skin
   # b) Buckling Skin (Mise à jour selon formule exacte)
    spacing = w_box / (n_str + 1) # Ceci correspond à 'b' sur l'image
    Kc = 4.0                      # Coefficient Kc pour bords simplement supportés
    
    # Formule exacte de l'image : Kc * [ (E * pi^2) / (12 * (1-v^2)) ] * (t/b)^2
    term_pre = (E_MODULUS * np.pi**2) / (12 * (1 - POISSON_RATIO**2))
    sigma_cr_skin = Kc * term_pre * (t_skin / spacing)**2
    
    margin_skin = (sigma_cr_skin / (SAFETY_FACTOR * sigma_f)) - 1.0
    
    # c) Buckling Stringer (Euler)
    L_rib = params['rib_spacing']
    I_str_own = (w_str * h_str**3) / 12 
    F_cr_str = (np.pi**2 * E_MODULUS * I_str_own) / (L_rib**2)
    sigma_cr_str = F_cr_str / (h_str * w_str)
    margin_str = (sigma_cr_str / (SAFETY_FACTOR * sigma_f)) - 1.0
    
    return [margin_yield, margin_skin, margin_str]

def wing_mass_objective_scaled(x_scaled, params):
    t_mm, h_cm, w_cm, n_str = x_scaled
    t_skin = t_mm / 1000.0
    h_str  = h_cm / 100.0
    w_str  = w_cm / 100.0
    
    b = params['span']
    
    Area_sect = (params['chord_root']*0.5 * t_skin) + (n_str * h_str * w_str)
    mass_wing_box = 2 * Area_sect * b * DENSITY * 0.7 
    return mass_wing_box * 1.25

def optimize_wing_structure(MTOW_current, geo_data):
    M_root = calculate_bending_moment_root(MTOW_current, geo_data['b'], geo_data['nz'])
    
    params = {
        'Moment_Flexion': M_root,
        'chord_root': geo_data['chord_root'],
        'tc_ratio': geo_data['tc'],
        'rib_spacing': 0.6,    
        'span': geo_data['b']
    }
    
    # Variables : [t_skin (mm), h_str (cm), w_str (cm), n_str (float)]
    # Point de départ un peu plus robuste
    x0 = [5.0, 10.0, 4.0, 30.0] 
    
    bounds = [
        (1.0, 60.0),    # Peau: jusqu'à 60mm
        (2.0, 40.0),    # Hauteur lisse
        (1.0, 15.0),    # Largeur lisse
        (2.0, 150.0)    # Nombre lisses (float bounds)
    ]
    
    cons = {'type': 'ineq', 'fun': lambda x: check_structure_constraints_scaled(x, params)}
    
    res = minimize(wing_mass_objective_scaled, x0, args=(params,), 
                   method='SLSQP', bounds=bounds, constraints=cons, 
                   options={'ftol': 1e-4, 'disp': False, 'maxiter': 100})
    
    # Conversion finale pour affichage
    x_final = res.x
    design_si = [x_final[0]/1000, x_final[1]/100, x_final[2]/100, x_final[3]] # keep n_str float for debug
    
    if res.success:
        return res.fun, design_si, True
    else:
        # On renvoie quand même le point final atteint pour analyse
        return res.fun, design_si, False

# 4. MAIN
def input_float(prompt, default_val):
    try:
        val_str = input(f"{prompt} [Défaut: {default_val}]: ")
        if val_str.strip() == "": return float(default_val)
        return float(val_str)
    except ValueError:
        return float(default_val)

if __name__ == "__main__":
    os.system('cls' if os.name == 'nt' else 'clear')
    print("="*70)
    print("   PROJET STRUCTURE AVION - BILAN DE MASSE & OPTIMISATION")
    print("="*70)
    
    # 1. ENTRÉES
    print("\n--- 1. CHARGEMENT ET GÉOMÉTRIE ---")
    mtow_in = input_float(" > MTOW cible [kg]", 77000.0)
    mzfw_in = input_float(" > MZFW (Zero Fuel Weight) [kg]", mtow_in * 0.78)
    b_in = input_float(" > Envergure (b) [m]", 34.1)
    c_root_in = input_float(" > Corde emplanture [m]", 6.0)
    tc_percent = input_float(" > Épaisseur relative [%]", 15.0)
    tc_in = tc_percent / 100.0
    taper_in = input_float(" > Effilement (Taper ratio) [0-1]", 0.3)
    sweep_in = input_float(" > Flèche à 25% (deg)", 30.0)
    
    S_calc = (b_in / 2) * c_root_in * (1 + taper_in)
    AR_calc = (b_in**2) / S_calc
    print(f"\n[GEOM] Surface: {S_calc:.1f} m² | Allongement: {AR_calc:.1f}")

    # 2. BILAN DE MASSE
    print("\n" + "-"*30)
    print("   BILAN DE MASSE (FORMULES THÈSE & RÉGRESSION)")
    print("-" * 30)

    m_kroo = formule_these_kroo(mtow_in, mzfw_in, S_calc, b_in, tc_in, taper_in, sweep_in)
    m_torenbeek = formule_these_torenbeek(mtow_in, mzfw_in, S_calc, b_in, tc_in, taper_in, sweep_in)
    res_reg = calcul_regression_custom(mtow_in)
    
    print(f"\n1. Thèse - Modèle I. Kroo (2001)     : {m_kroo:.0f} kg")
    print(f"2. Thèse - Modèle E. Torenbeek (1986): {m_torenbeek:.0f} kg")
    
    if res_reg:
        m_custom, a, b_exp = res_reg
        print(f"3. Votre Création (Régression)       : {m_custom:.0f} kg  (Formule: {a:.4f} * MTOW^{b_exp:.4f})")
    else:
        m_custom = None
        print("3. Votre Création                    : INDISPONIBLE")

    # 3. OPTIMISATION
    print("\n" + "-"*30)
    print("   OPTIMISATION STRUCTURELLE (CALCUL PHYSIQUE)")
    print("-" * 30)
    
    geo_data = {'b': b_in, 'chord_root': c_root_in, 'tc': tc_in, 'nz': 3.75}
    
    print("... Optimisation des sections (Peau/Lisses) en cours ...")
    mass_opt, design_opt, success = optimize_wing_structure(mtow_in, geo_data)
    
    # Analyse marges finales
    M_final = calculate_bending_moment_root(mtow_in, b_in, 3.75)
    params_check = {'Moment_Flexion': M_final, 'chord_root': c_root_in, 'tc_ratio': tc_in, 'rib_spacing': 0.6, 'span': b_in}
    
    # Rescaling pour le check
    # design_opt est en SI [m, m, m, float]
    # check attend [mm, cm, cm, float]
    x_check = [design_opt[0]*1000, design_opt[1]*100, design_opt[2]*100, design_opt[3]]
    margs = check_structure_constraints_scaled(x_check, params_check)

    print("-" * 60)
    if success:
        print(f"RÉSULTAT SUCCÈS : MASSE AILE CALCULÉE = {mass_opt:.1f} kg")
    else:
        print(f"/!\\ ATTENTION : L'optimiseur n'a pas convergé parfaitement.")
        print(f"Meilleur point trouvé : {mass_opt:.1f} kg")

    # Comparaisons
    diff_kroo = (mass_opt - m_kroo) / m_kroo * 100
    print(f" > Écart vs Kroo      : {diff_kroo:+.1f} %")
    diff_tor = (mass_opt - m_torenbeek) / m_torenbeek * 100
    print(f" > Écart vs Torenbeek : {diff_tor:+.1f} %")
    if m_custom:
        diff_reg = (mass_opt - m_custom) / m_custom * 100
        print(f" > Écart vs Régression: {diff_reg:+.1f} %")

    t_mm = design_opt[0] * 1000
    h_mm = design_opt[1] * 1000
    w_mm = design_opt[2] * 1000
    n_str = int(design_opt[3]) # Affichage en entier
    
    print("\n>>> Configuration retenue (Niveau 1 & 2):")
    print(f"  Peau (t)     : {t_mm:.2f} mm")
    print(f"  Lisses       : {n_str} lisses de section {h_mm:.1f}x{w_mm:.1f} mm")
    
    print("\n>>> Marges de sécurité (Objectif >= 0.0):")
    print(f"  Plasticité (Yield) : {margs[0]:+.2f}  {'[OK]' if margs[0]>= -0.01 else '[ECHEC]'}")
    print(f"  Flambage Peau      : {margs[1]:+.2f}  {'[OK]' if margs[1]>= -0.01 else '[ECHEC]'}")
    print(f"  Flambage Lisses    : {margs[2]:+.2f}  {'[OK]' if margs[2]>= -0.01 else '[ECHEC]'}")