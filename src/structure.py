import numpy as np
import os
import sys

# Import du module de régression (assurez-vous que regression.py est dans le même dossier)
try:
    from regression import effectuer_regression_puissance
except ImportError:
    print("ATTENTION : Le fichier regression.py est introuvable.")
    # Valeurs par défaut si le fichier manque (backup)
    effectuer_regression_puissance = lambda afficher_graphe: (0.0269, 1.1223)

#CONSTANTES MATÉRIAUX (Aluminium 2024-T3 / 7075-T6)
E_ALU = 73e9       # Module d'Young (Pa)
RHO_ALU = 2780     # Masse volumique (kg/m3)
SIGMA_Y = 450e6    # Limite élastique (Pa) - Yield Stress
NU_ALU = 0.33      # Coefficient de Poisson

class AircraftStructure:
    """
    Outil de dimensionnement structurel pour le projet VAD.
    Permet d'estimer les masses (Ailes, Fuselage) et de vérifier la tenue mécanique.
    """

    def __init__(self, factor_safety=1.5):
        """
        Args:
            factor_safety (float): Coefficient de sécurité (1.5 standard aéro).
        """
        self.sf = factor_safety
        
        # RÉCUPÉRATION DYNAMIQUE DES COEFFICIENTS DE RÉGRESSION
        print("Chargement du modèle de régression...")
        # On met afficher_graphe=False pour ne pas bloquer le script avec une fenêtre
        self.reg_a, self.reg_b = effectuer_regression_puissance(afficher_graphe=False)
        
        print(f"Modèle chargé : M_wing = {self.reg_a:.4f} * MTOW ^ {self.reg_b:.4f}")

    # 1. ESTIMATION DE MASSE : VOILURE
    def wing_mass_regression_computed(self, MTOW_kg):
        """
        ÉQUATION DE CRÉATION (Groupe 1) :
        Utilise la Loi de Puissance calculée dynamiquement par regression.py.
        Formule : M_wing = a * MTOW ^ b
        """
        if self.reg_a == 0 or self.reg_b == 0:
            print("AVERTISSEMENT : Régression invalide (coeffs nuls).")
            return 0.0
            
        return self.reg_a * (MTOW_kg ** self.reg_b)

    def wing_mass_roux_secondary(self, S_m2, is_jet=True, mtow_tonnes=70):
        """
        FORMULE 1 (Thèse E. Roux, p.167) : Structure Secondaire
        (Becs, volets, aérofreins, carénages...)
        """
        k_helice = 1.0 if is_jet else 0.488
        
        if mtow_tonnes >= 20:
            K = 25.9
            n = 0.97
        else:
            K = 4.39
            n = 1.358
            
        M_wss = k_helice * K * (S_m2 ** n)
        return M_wss

    def wing_mass_roux_other(self, M_caisson_kg):
        """
        FORMULE 2 (Thèse E. Roux, p.175) : Autres éléments
        """
        return 0.12 * M_caisson_kg

    def wing_mass_leclerc_primary(self, MTOW_kg, S_m2, b_m, phi_rad, t_c=0.12):
        """
        Formule de F. Leclerc simplifiée pour estimer la masse du CAISSON.
        """
        AR = b_m**2 / S_m2 # Allongement
        epsilon = 0.3
        
        term = (MTOW_kg * AR * (epsilon + 1)) / (t_c * np.cos(phi_rad))
        M_caisson = 0.15 * (term**0.6) * (S_m2**0.3)
        return M_caisson

    def compute_wing_mass_analytic(self, MTOW_kg, S_m2, b_m, phi_rad):
        """
        Méthode Analytique Complète (Hybride Leclerc + Roux).
        """
        # 1. Masse Caisson (Leclerc)
        m_prim = self.wing_mass_leclerc_primary(MTOW_kg, S_m2, b_m, phi_rad)
        
        # 2. Masse Secondaire (Roux Formule 1)
        m_sec = self.wing_mass_roux_secondary(S_m2, is_jet=True, mtow_tonnes=MTOW_kg/1000)
        
        # 3. Masse Autres (Roux Formule 2)
        m_autre = self.wing_mass_roux_other(m_prim)
        
        return m_prim + m_sec + m_autre

    # 2. ESTIMATION DE MASSE : FUSELAGE
    def fuselage_mass_kroo(self, L_fus_m, D_fus_m, MTOW_kg, delta_P_Pa=55000):
        """
        Modèle d'Ilan Kroo (Stanford).
        """
        L_ft = L_fus_m * 3.28084
        D_ft = D_fus_m * 3.28084
        MTOW_lb = MTOW_kg * 2.20462
        dP_lbft2 = delta_P_Pa * 0.020885 

        Ip = 1.5e-3 * dP_lbft2 * D_ft
        
        nz = 2.5 
        m_ch = 0.53 * MTOW_lb 
        Ib = 1.91e-4 * nz * m_ch * (L_ft / (D_ft**2))
        
        if Ip > Ib:
            If = Ip
        else:
            If = (Ip**2 + Ib**2) / (2*Ib)
            
        Sw_ft2 = np.pi * D_ft * L_ft * ((1 - 2*D_ft/L_ft)**(2/3)) * (1 + (D_ft/L_ft)**2)
        
        M_fus_lb = (1.051 + 0.102 * If) * Sw_ft2
        
        return M_fus_lb / 2.20462 

    # 3. VÉRIFICATION MÉCANIQUE (RDM)
    def check_fuselage_integrity(self, R_fus_m, t_skin_m, delta_P_Pa, Moment_Flexion_Nm):
        """
        Vérifie si l'épaisseur du fuselage (t_skin) est suffisante.
        """
        # 1. Contrainte de Pression
        sigma_hoop = (delta_P_Pa * R_fus_m) / t_skin_m
        
        # 2. Contrainte de Flexion
        I_yy = np.pi * (R_fus_m**3) * t_skin_m
        sigma_bend = (Moment_Flexion_Nm * R_fus_m) / I_yy
        
        # 3. Critère de Rupture
        sigma_total = sigma_hoop + sigma_bend
        
        # Marge de sécurité
        margin = SIGMA_Y - (sigma_total * self.sf)
        is_safe = margin > 0
        
        return {
            "safe": is_safe,
            "margin_MPa": margin / 1e6,
            "sigma_hoop_MPa": sigma_hoop / 1e6,
            "sigma_bend_MPa": sigma_bend / 1e6,
            "sigma_total_MPa": sigma_total / 1e6
        }

    def check_buckling_spar(self, L_element, F_comp, h, b, e):
        """
        Vérification flambage longeron.
        """
        I = (e * h**3) / 12 
        F_crit = (np.pi**2 * E_ALU * I) / (L_element**2)
        return F_crit > (F_comp * self.sf), F_crit


# TEST DU MODULE
if __name__ == "__main__":
    print("\n--- TEST DU MODULE STRUCTURE (PROJET VAD) ---\n")
    
    # 1. Initialisation (Cela va lancer la régression en arrière-plan)
    tool = AircraftStructure()
    
    # 2. Définition d'un Avion "Hypothèse"
    hypothese_MTOW = 72000.0  # kg
    S_aile = 122.6            # m2
    b_aile = 34.1             # m
    phi_aile = np.radians(25) 
    
    L_fus = 37.6              # m
    D_fus = 3.95              # m
    R_fus = D_fus / 2
    
    # 3. Calculs de Masse
    print("[1] ESTIMATION DES MASSES")
    
    # Méthode Régression (DYNAMIQUE MAINTENANT)
    m_wing_reg = tool.wing_mass_regression_computed(hypothese_MTOW)
    print(f"  - Masse Aile (Loi Puissance Calculée) : {m_wing_reg:.1f} kg")
    
    # Méthode Analytique (Roux + Leclerc)
    m_wing_ana = tool.compute_wing_mass_analytic(hypothese_MTOW, S_aile, b_aile, phi_aile)
    print(f"  - Masse Aile (Analytique Roux)        : {m_wing_ana:.1f} kg")
    
    # Masse Fuselage (Kroo)
    m_fus = tool.fuselage_mass_kroo(L_fus, D_fus, hypothese_MTOW)
    print(f"  - Masse Fuselage (Kroo)               : {m_fus:.1f} kg")
    
    print("-" * 40)
    
    # 4. Vérification Structurelle (Fuselage)
    print("[2] VÉRIFICATION FUSELAGE")
    
    delta_P = 55000.0       
    M_flexion = 2.5e6       
    t_peau = 0.0018         
    
    res = tool.check_fuselage_integrity(R_fus, t_peau, delta_P, M_flexion)
    
    print(f"  - Épaisseur testée : {t_peau*1000} mm")
    print(f"  - Contrainte Pression : {res['sigma_hoop_MPa']:.1f} MPa")
    print(f"  - Contrainte Flexion  : {res['sigma_bend_MPa']:.1f} MPa")
    print(f"  - Total (Von Mises)   : {res['sigma_total_MPa']:.1f} MPa")
    
    if res['safe']:
        print(f"  => RÉSULTAT : CONFORME (Marge : {res['margin_MPa']:.1f} MPa)")
    else:
        print(f"  => RÉSULTAT : ÉCHEC (Rupture ou Plastification)")