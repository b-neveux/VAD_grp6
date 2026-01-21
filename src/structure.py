import numpy as np
import os
import sys

# Import du module de régression
try:
    from regression import effectuer_regression_puissance
except ImportError:
    print("ATTENTION : Le fichier regression.py est introuvable.")
    effectuer_regression_puissance = lambda afficher_graphe: (0.0269, 1.1223)

# CONSTANTES MATÉRIAUX
E_ALU = 73e9       
RHO_ALU = 2780    
SIGMA_Y = 450e6    
NU_ALU = 0.33      

class AircraftStructure:
    """
    Outil de dimensionnement structurel pour le projet VAD.
    """
    def __init__(self, factor_safety=1.5):
        self.sf = factor_safety
        print("\n--- CHARGEMENT ---")
        print("Chargement du modèle de régression statistique...")
        self.reg_a, self.reg_b = effectuer_regression_puissance(afficher_graphe=False)
        # On supprime le print ici pour garder la console propre pour les questions

    # 1. ESTIMATION DE MASSE : VOILURE
    def wing_mass_leclerc_primary(self, MTOW_kg, S_m2, b_m, phi_rad, t_c):
        """
        Note: t_c doit être passé explicitement ici (ex: 0.15 pour 15%)
        """
        AR = b_m**2 / S_m2 
        epsilon = 0.3
        # Protection contre la division par zéro si t_c est mal rentré
        if t_c < 0.01: t_c = 0.12 
        
        term = (MTOW_kg * AR * (epsilon + 1)) / (t_c * np.cos(phi_rad))
        M_caisson = 0.15 * (term**0.6) * (S_m2**0.3)
        return M_caisson

    def wing_mass_roux_secondary(self, S_m2, is_jet=True, mtow_tonnes=70):
        k_helice = 1.0 if is_jet else 0.488
        if mtow_tonnes >= 20:
            K = 25.9; n = 0.97
        else:
            K = 4.39; n = 1.358
        return k_helice * K * (S_m2 ** n)

    def wing_mass_roux_other(self, M_caisson_kg):
        return 0.12 * M_caisson_kg

    def compute_wing_mass_analytic(self, MTOW_kg, S_m2, b_m, phi_rad, t_c_ratio):
        # On passe t_c_ratio (épaisseur) à la formule de Leclerc
        m_prim = self.wing_mass_leclerc_primary(MTOW_kg, S_m2, b_m, phi_rad, t_c=t_c_ratio)
        m_sec = self.wing_mass_roux_secondary(S_m2, is_jet=True, mtow_tonnes=MTOW_kg/1000)
        m_autre = self.wing_mass_roux_other(m_prim)
        return m_prim + m_sec + m_autre

    # 2. ESTIMATION DE MASSE : FUSELAGE
    def fuselage_mass_kroo(self, L_fus_m, D_fus_m, MTOW_kg, delta_P_Pa=55000):
        L_ft = L_fus_m * 3.28084
        D_ft = D_fus_m * 3.28084
        MTOW_lb = MTOW_kg * 2.20462
        dP_lbft2 = delta_P_Pa * 0.020885 

        Ip = 1.5e-3 * dP_lbft2 * D_ft
        nz = 2.5 
        m_ch = 0.53 * MTOW_lb 
        Ib = 1.91e-4 * nz * m_ch * (L_ft / (D_ft**2))
        
        If = Ip if Ip > Ib else (Ip**2 + Ib**2) / (2*Ib)
        Sw_ft2 = np.pi * D_ft * L_ft * ((1 - 2*D_ft/L_ft)**(2/3)) * (1 + (D_ft/L_ft)**2)
        M_fus_lb = (1.051 + 0.102 * If) * Sw_ft2
        return M_fus_lb / 2.20462 

    # 3. RDM
    def check_fuselage_integrity(self, R_fus_m, t_skin_m, delta_P_Pa, Moment_Flexion_Nm):
        sigma_hoop = (delta_P_Pa * R_fus_m) / t_skin_m
        I_yy = np.pi * (R_fus_m**3) * t_skin_m
        sigma_bend = (Moment_Flexion_Nm * R_fus_m) / I_yy
        sigma_total = sigma_hoop + sigma_bend
        margin = SIGMA_Y - (sigma_total * self.sf)
        return {"safe": margin > 0, "margin_MPa": margin / 1e6}

    # 4. EXPORT SIMULATEUR
    def estimer_inerties_et_export(self, MTOW_kg, b_aile, L_fus, m_aile, m_fus):
        # Répartition de masse simplifiée
        m_reste = MTOW_kg - m_aile - m_fus
        
        x_fus = 0.45 * L_fus      
        x_aile = 0.40 * L_fus     
        x_reste = 0.45 * L_fus    
        
        x_cg = (m_fus * x_fus + m_aile * x_aile + m_reste * x_reste) / MTOW_kg
        
        # Rayons de giration approximatifs (Raymer)
        Rx = 0.25 * b_aile
        Ixx = MTOW_kg * (Rx ** 2)
        
        Ry = 0.34 * L_fus
        Iyy = MTOW_kg * (Ry ** 2)
        
        Izz = Ixx + Iyy # Approx simple
        
        return {
            "m": round(MTOW_kg, 2),
            "Ixx": round(Ixx, 2), "Iyy": round(Iyy, 2), "Izz": round(Izz, 2),
            "Ixy": 0, "Ixz": 0, "Iyz": 0,
            "CG_pos_m": round(x_cg, 2)
        }

# --- FONCTION UTILITAIRE POUR DEMANDER À L'UTILISATEUR ---
def input_float(prompt, default_val):
    try:
        val_str = input(f"{prompt} [Défaut: {default_val}]: ")
        if val_str.strip() == "":
            return float(default_val)
        return float(val_str)
    except ValueError:
        print(f"Erreur de saisie. Utilisation de la valeur par défaut : {default_val}")
        return float(default_val)

# --- MAIN INTERACTIF ---
if __name__ == "__main__":
    os.system('cls' if os.name == 'nt' else 'clear') # Nettoyer la console
    print("\n" + "="*60)
    print("   EXTRACTION XFLR5 -> SIMULATEUR DE VOL")
    print("   Entrez les valeurs lues dans XFLR5 (Current Plane -> Define)")
    print("="*60 + "\n")
    
    tool = AircraftStructure()
    
    # 1. SAISIE DES DONNÉES UTILISATEUR
    print("\n--- 1. GÉOMÉTRIE & MASSE CIBLE ---")
    mtow_in = input_float("Masse Totale au décollage (MTOW) en kg ?", 72000)
    s_aile_in = input_float("Surface Alaire (S) en m² ?", 122.6)
    b_aile_in = input_float("Envergure (b) en m ?", 34.1)
    mac_in = input_float("Corde Moyenne (MAC) donnée par XFLR5 en m ?", 3.8) # Mieux que S/b
    
    print("\n--- 2. PARAMÈTRES AVANCÉS (Pour la structure) ---")
    print("Info : Regardez votre profil dans XFLR5 (ex: NACA 2415 -> 15%)")
    tc_percent = input_float("Épaisseur relative du profil (%) ?", 12.0)
    tc_val = tc_percent / 100.0 # Conversion en décimal (0.12)
    
    phi_deg = input_float("Angle de flèche de l'aile (degrés) ?", 25.0)
    phi_rad = np.radians(phi_deg)
    
    l_fus_in = input_float("Longueur du fuselage (m) ?", 37.6)
    d_fus_in = input_float("Diamètre du fuselage (m) ?", 3.95)

    # 2. CALCULS
    print("\n" + "-"*30)
    print("Traitement en cours...")
    
    # On passe le tc_val (épaisseur) au calcul de masse !
    m_wing = tool.compute_wing_mass_analytic(mtow_in, s_aile_in, b_aile_in, phi_rad, tc_val)
    m_fus = tool.fuselage_mass_kroo(l_fus_in, d_fus_in, mtow_in)
    
    sim_data = tool.estimer_inerties_et_export(mtow_in, b_aile_in, l_fus_in, m_wing, m_fus)
    
    # 3. AFFICHAGE FINAL POUR LE SIMULATEUR
    print("\n" + "="*60)
    print("   RÉSULTATS À COPIER DANS LE SIMULATEUR")
    print("   (Onglet 1 : Géométrie, masses et inerties)")
    print("="*60)
    
    # Affichage aligné pour lecture facile
    print(f"{'VARIABLE':<10} | {'VALEUR':<15} | {'UNITÉ':<5}")
    print("-" * 35)
    print(f"{'c':<10} | {mac_in:<15.3f} | m")
    print(f"{'b':<10} | {b_aile_in:<15.3f} | m")
    print(f"{'S':<10} | {s_aile_in:<15.3f} | m²")
    print(f"{'m':<10} | {sim_data['m']:<15.1f} | kg")
    print("-" * 35)
    print(f"{'Ixx':<10} | {sim_data['Ixx']:<15.0f} | kg.m²")
    print(f"{'Iyy':<10} | {sim_data['Iyy']:<15.0f} | kg.m²")
    print(f"{'Izz':<10} | {sim_data['Izz']:<15.0f} | kg.m²")
    print(f"{'Ixy':<10} | {sim_data['Ixy']:<15.0f} | -")
    print(f"{'Ixz':<10} | {sim_data['Ixz']:<15.0f} | -")
    print(f"{'Iyz':<10} | {sim_data['Iyz']:<15.0f} | -")
    print("="*60)
    
    # Petit bilan de masse structurelle pour info
    print(f"\n[INFO STRUCTURE] Masse Aile estimée : {m_wing:.0f} kg (Profil {tc_percent}%)")
    print(f"[INFO STRUCTURE] Masse Fuselage est.  : {m_fus:.0f} kg")