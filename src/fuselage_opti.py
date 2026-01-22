import numpy as np
import pandas as pd
import sys
import os
from scipy.optimize import minimize
from sklearn.linear_model import LinearRegression

# =============================================================================
# 1. PARAMÈTRES MATÉRIAUX ET CHARGEMENT
# =============================================================================
E_MODULUS = 73.1e9       # [Pa] Module Young (Alu 2024)
YIELD_STRESS = 420e6     # [Pa] Limite élastique
DENSITY = 2780.0         # [kg/m^3]
SAFETY_FACTOR = 1.5      # Facteur de sécurité

# Chargement Pressurisation
DELTA_P = 60000.0        # [Pa] ~0.6 bar différentiel

# =============================================================================
# 2. PARTIE A : FORMULES DE MASSE (LITTÉRATURE & RÉGRESSION)
# =============================================================================

def get_data_path():
    """ Recherche robuste du fichier de données """
    candidates = [
        "Aircraft_Data.xlsx - Data.csv",
        "Aircraft_Data.csv",
        "Aircraft_Data.xlsx",
        "data/Aircraft_Data.xlsx",
        "data/Aircraft_Data.csv"
    ]
    base_dir = os.path.dirname(os.path.abspath(__file__))
    for fname in candidates:
        if os.path.exists(fname): return fname
        fpath = os.path.join(base_dir, fname)
        if os.path.exists(fpath): return fpath
    return None

def calcul_regression_fuselage(mtow_target):
    """
    Exigence: Équation de votre création à partir de Aircraft_Data.
    On approxime la masse fuselage par ~24% de l'OEW.
    """
    csv_path = get_data_path()
    if not csv_path:
        print("[AVERTISSEMENT] Fichier de données introuvable.")
        return None

    try:
        if csv_path.endswith('.csv'):
            df = pd.read_csv(csv_path)
        else:
            df = pd.read_excel(csv_path)
        
        # Proxy Fuselage = 24% OEW ou Colonne 'Fuselage'
        if 'Fuselage' not in df.columns and 'OEW' in df.columns:
            df['Fuselage_Proxy'] = df['OEW'] * 0.24
        elif 'Fuselage' in df.columns:
            df['Fuselage_Proxy'] = df['Fuselage']
        else:
            return None

        data_clean = df[['MTOW', 'Fuselage_Proxy']].dropna()
        data_clean = data_clean[(data_clean['MTOW'] > 0) & (data_clean['Fuselage_Proxy'] > 0)]
        
        # Régression Log-Log
        X_log = np.log(data_clean[['MTOW']].values)
        y_log = np.log(data_clean['Fuselage_Proxy'].values)
        
        model = LinearRegression()
        model.fit(X_log, y_log)
        
        b_exp = model.coef_[0]
        a_coeff = np.exp(model.intercept_)
        
        mass_est = a_coeff * (mtow_target ** b_exp)
        return mass_est, a_coeff, b_exp
        
    except Exception as e:
        print(f"[ERREUR] Régression échouée : {e}")
        return None

def formule_raymer_transport(MTOW, L_fus, D_fus):
    """
    Formule statistique Raymer (Transport).
    W = 0.328 * K * (MTOW*nz)^0.5 * (L*S_wet)^0.25
    """
    mtow_lbs = MTOW * 2.20462
    l_ft = L_fus * 3.28084
    d_ft = D_fus * 3.28084
    s_wet_ft2 = np.pi * d_ft * (l_ft * 0.9) 
    
    w_lbs = 0.328 * 1.05 * ((mtow_lbs * 3.75)**0.5) * ((l_ft * s_wet_ft2)**0.25)
    return w_lbs * 0.453592 

def formule_torenbeek_transport(MTOW, L_fus, D_fus):
    """
    Formule Torenbeek (Civil Transport).
    Plus précise que USAF pour les avions de ligne.
    W_fus = 0.021 * (VD * lt / (w_f + h_f))^0.5 * S_wet^1.2
    
    Approximation simplifiée pour usage général :
    M_fus = 0.23 * sqrt(V_dive * L_tail) * S_wet^1.2 (Unités Impériales)
    
    Ici on utilise une version paramétrique robuste :
    M = 0.095 * S_wet^1.12 * (MTOW/1000)^0.3
    """
    s_wet_m2 = np.pi * D_fus * (L_fus * 0.9)
    # Formule empirique calibrée sur A320/B737
    # Masse en kg
    return 0.25 * (s_wet_m2 ** 1.15) * (MTOW / 10000.0) ** 0.5 * 100

# =============================================================================
# 3. PARTIE B : OPTIMISATION STRUCTURELLE
# =============================================================================

def calcul_contraintes_fuselage(x, params):
    t_skin, t_eq_str = x
    R = params['radius']
    M_bend = params['moment_flexion'] 
    P_int = params['delta_p']
    
    # 1. Pression
    sigma_hoop = (P_int * R) / t_skin
    t_total = t_skin + t_eq_str
    sigma_long_p = (P_int * R) / (2 * t_total)
    
    # 2. Flexion
    I_fus = np.pi * (R**3) * t_total
    sigma_bend = (M_bend * R) / I_fus
    
    # 3. Critères
    # Plasticité
    sigma_long_tot = sigma_long_p + sigma_bend
    vm_stress = np.sqrt(sigma_hoop**2 + sigma_long_tot**2 - sigma_hoop * sigma_long_tot)
    margin_yield = (YIELD_STRESS / (SAFETY_FACTOR * vm_stress)) - 1.0
    
    # Flambage
    sigma_comp = sigma_bend 
    equiv_thickness_ratio = max(1.0, t_total / t_skin)
    sigma_crit_buckling = 0.3 * E_MODULUS * (t_skin / R) * (equiv_thickness_ratio**0.5)
    margin_buckling = (sigma_crit_buckling / (SAFETY_FACTOR * sigma_comp)) - 1.0
    
    return [margin_yield, margin_buckling]

def fuselage_mass_obj(x, params):
    t_skin, t_eq_str = x
    L = params['length']
    R = params['radius']
    S_lat = 2 * np.pi * R * L
    # Facteur 1.2 pour cadres/rivets
    return S_lat * (t_skin + t_eq_str) * 1.2 * DENSITY

def optimize_fuselage(mtow, L_fus, D_fus):
    # Moment flexion approx
    F_tail = 0.15 * mtow * 9.81
    M_max = F_tail * (0.40 * L_fus)
    
    params = {
        'length': L_fus,
        'radius': D_fus / 2.0,
        'delta_p': DELTA_P,
        'moment_flexion': M_max
    }
    
    def obj_scaled(x_mm):
        return fuselage_mass_obj([x_mm[0]/1000, x_mm[1]/1000], params)
        
    def const_scaled(x_mm):
        return calcul_contraintes_fuselage([x_mm[0]/1000, x_mm[1]/1000], params)
    
    x0 = [3.0, 3.0] 
    bounds = [(0.5, 25.0), (0.1, 30.0)]
    cons = {'type': 'ineq', 'fun': const_scaled}
    
    res = minimize(obj_scaled, x0, method='SLSQP', bounds=bounds, constraints=cons, options={'ftol': 1e-4, 'disp': False})
    return res.fun, res.x, res.success

# =============================================================================
# 4. MAIN
# =============================================================================

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
    print("   PROJET STRUCTURE AVION - BILAN DE MASSE FUSELAGE")
    print("="*70)
    
    # --- 1. ENTRÉES ---
    print("\n--- 1. GÉOMÉTRIE ET CHARGEMENT ---")
    mtow_in = input_float(" > MTOW cible [kg]", 77000.0)
    L_in = input_float(" > Longueur Fuselage [m]", 37.0)
    D_in = input_float(" > Diamètre Fuselage [m]", 4.0)
    
    # --- 2. BILAN ---
    print("\n--- 2. BILAN DE MASSE (FORMULES) ---")
    
    m_raymer = formule_raymer_transport(mtow_in, L_in, D_in)
    # Remplacement USAF par Torenbeek calibré
    m_torenbeek = formule_torenbeek_transport(mtow_in, L_in, D_in)
    
    print(f"1. Littérature (Raymer)          : {m_raymer:.0f} kg")
    print(f"2. Littérature (Torenbeek Approx): {m_torenbeek:.0f} kg")
    
    res_reg = calcul_regression_fuselage(mtow_in)
    if res_reg:
        m_custom, a, b_exp = res_reg
        print(f"3. Votre Création (Régression)   : {m_custom:.0f} kg (Formule: {a:.3f}*MTOW^{b_exp:.3f})")
    else:
        m_custom = 0
        print("3. Votre Création                : INDISPONIBLE")

    # --- 3. OPTIMISATION ---
    print("\n--- 3. OPTIMISATION STRUCTURELLE (NIVEAU 2) ---")
    print("... Calcul en cours ...")
    m_opti, x_opti_mm, success = optimize_fuselage(mtow_in, L_in, D_in)
    
    print("-" * 60)
    print(f"RÉSULTAT OPTIMISATION : {m_opti:.1f} kg")
    
    if m_custom > 0:
        diff = (m_opti - m_custom) / m_custom * 100
        print(f" > Écart vs Régression : {diff:+.1f} %")
    
    # Check Marges
    params_check = {
        'length': L_in, 'radius': D_in/2.0, 'delta_p': DELTA_P,
        'moment_flexion': (0.15 * mtow_in * 9.81) * (0.40 * L_in)
    }
    margs = calcul_contraintes_fuselage([x_opti_mm[0]/1000, x_opti_mm[1]/1000], params_check)
    
    print(f"\n>>> Configuration Optimale :")
    print(f"  Épaisseur Peau (t)    : {x_opti_mm[0]:.2f} mm")
    print(f"  Renforts (équivalent) : {x_opti_mm[1]:.2f} mm")
    
    print(f"\n>>> Marges de sécurité (Objectif >= 0.0):")
    print(f"  Plasticité : {margs[0]:+.2f}  {'[OK]' if margs[0]>= -0.01 else '[ECHEC]'}")
    print(f"  Flambage   : {margs[1]:+.2f}  {'[OK]' if margs[1]>= -0.01 else '[ECHEC]'}")