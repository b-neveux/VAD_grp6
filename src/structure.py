import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import optimization modules
try:
    import wing_opti
    import fuselage_opti
except ImportError:
    print("[ERREUR] Les fichiers 'wing_opti.py' et 'fuselage_opti.py' doivent être dans le même dossier.")
    sys.exit(1)

# =============================================================================
# CONSTANTES STATISTIQUES (Ratios de masse Roskam/Torenbeek)
# =============================================================================
RATIO_EMPENNAGE = 0.025   # ~2.5% MTOW
RATIO_TRAIN     = 0.040   # ~4.0% MTOW
RATIO_NACELLES  = 0.018   # ~1.8% MTOW
RATIO_SYSTEMES  = 0.120   # ~12% MTOW (Hydraulique, Elec, Avionique...)
RATIO_AMENAGEMENT= 0.060  # ~6% MTOW (Sièges, Galleys...)

# Constantes Physiques & Mission
G = 9.81
V_SOUND_35K = 295.0       # Vitesse son à 35000ft [m/s]
DIST_LILLE_COPENHAGUE = 800.0 # [km] Distance opérationnelle estimée

def input_float(prompt, default_val):
    try:
        val_str = input(f"{prompt} [Défaut: {default_val}]: ")
        if val_str.strip() == "": return float(default_val)
        return float(val_str)
    except ValueError:
        return float(default_val)

def calcul_carburant_breguet(mtow, range_km, mach, sfc_mg, lod):
    """
    Calcule la masse de carburant nécessaire via l'équation de Bréguet.
    R = (V / (g * SFC)) * (L/D) * ln(Wi / Wf)
    """
    range_m = range_km * 1000.0
    velocity = mach * V_SOUND_35K  # m/s
    sfc_si = sfc_mg * 1e-6         # mg/N/s -> kg/N/s
    
    # W_final / W_initial = exp( - (R * g * SFC) / (V * L/D) )
    exponent = (range_m * G * sfc_si) / (velocity * lod)
    mass_ratio = np.exp(-exponent) 
    
    fuel_fraction = 1.0 - mass_ratio
    mission_fuel = mtow * fuel_fraction
    total_fuel = mission_fuel * 1.05 # +5% Réserves
    
    return total_fuel, mission_fuel

# =============================================================================
# FONCTIONS D'INERTIE ET GEOMETRIE AVANCÉE
# =============================================================================

def calcul_inertie_globale(m_comps, geo):
    """
    Estime le tenseur d'inertie de l'avion complet.
    Référentiel Avion : X vers l'arrière, Y vers la droite, Z vers le bas.
    Origine : Nez de l'appareil (X=0, Y=0, Z=0).
    """
    # 1. Récupération des masses
    m_wing = m_comps['wing']
    m_fus = m_comps['fuselage']
    m_tail = m_comps['empennage']
    m_gear = m_comps['train']
    m_prop = m_comps['propulsion'] # Moteurs + Nacelles
    m_sys = m_comps['systemes']    # Systèmes + Aménagements
    m_pay = m_comps['payload']
    m_fuel = m_comps['fuel']

    # 2. Géométrie
    b = geo['b']
    c_root = geo['c_root']
    taper = geo['taper']
    sweep_25 = geo['sweep_25']
    L_fus = geo['L_fus']
    D_fus = geo['D_fus']

    # 3. Hypothèses de Position (Centres de Gravité Locaux)
    # Positions X (Longitudinal, depuis le nez)
    X_fus_cg = 0.45 * L_fus
    X_wing_apex = 0.35 * L_fus  # Bord d'attaque emplanture à 35% du fuselage
    X_tail = 0.94 * L_fus       # Empennage en queue
    X_eng = X_wing_apex + 0.5 * c_root # Moteurs sous l'aile (approx mi-corde)
    X_gear = X_wing_apex + 0.6 * c_root
    X_sys = 0.45 * L_fus        # Réparti
    X_pay = 0.45 * L_fus        # Réparti
    
    # Positions Z (Vertical, Z=0 à l'axe fuselage, >0 vers le bas)
    Z_fus = 0.0
    Z_wing = 0.2 * D_fus    # Aile basse
    Z_eng = 0.5 * D_fus     # Moteurs sous l'aile
    Z_tail = -0.3 * D_fus   # Empennage un peu surélevé
    Z_gear = 0.8 * D_fus    # Train sorti (ou rétracté ~0.2) -> Prenons moyenné 0.5
    Z_pay = 0.0
    
    # --- CALCUL INERTIE VOILURE (STRIP THEORY) ---
    # On intègre aile + fuel
    mass_wing_fuel = m_wing + m_fuel
    
    # Conversion flèche 25% -> Flèche Bord d'Attaque (LE)
    c_tip = c_root * taper
    tan_phi_25 = np.tan(np.radians(sweep_25))
    x_LE_tip = (b/2) * tan_phi_25 + 0.25 * (c_tip - c_root)
    sweep_LE_rad = np.arctan(x_LE_tip / (b/2))
    
    # Discrétisation
    N = 50
    dy = (b/2) / N
    y_strips = np.linspace(dy/2, b/2 - dy/2, N)
    
    # Listes pour sommer
    mass_strips = []
    
    total_area_strip = 0
    for y in y_strips:
        c_y = c_root * (1 - (1-taper) * (y / (b/2)))
        total_area_strip += c_y * dy

    # Calcul propriétés aile (par rapport au Nez)
    # I_wing_nose (Tensur 3x3 à l'origine Nez)
    Ixx_w, Iyy_w, Izz_w, Ixz_w = 0, 0, 0, 0
    
    MX_w, MY_w, MZ_w = 0, 0, 0 # Moments statiques pour trouver CG aile
    
    for y in y_strips:
        c_y = c_root * (1 - (1-taper) * (y / (b/2)))
        dm = (mass_wing_fuel / (total_area_strip * 2)) * (c_y * dy) * 2 # *2 pour symétrie G/D
        
        # Position locale
        x_le_local = y * np.tan(sweep_LE_rad)
        x_loc = X_wing_apex + x_le_local + 0.4 * c_y # CG structural à 40% corde
        z_loc = Z_wing
        # y_loc = y (et -y)
        
        # Contribution Inertie (Théorème Huygens pour point mass)
        Ixx_w += dm * (y**2 + z_loc**2)
        Iyy_w += dm * (x_loc**2 + z_loc**2)
        Izz_w += dm * (x_loc**2 + y**2)
        Ixz_w += dm * (x_loc * z_loc)
        
        MX_w += dm * x_loc
        MZ_w += dm * z_loc
        # MY_w = 0 (Symétrie)

    # --- AUTRES COMPOSANTS (POINTS MASSES OU CYLINDRES) ---
    components = [
        {'m': m_fus, 'x': X_fus_cg, 'y': 0, 'z': Z_fus, 'type': 'cyl_x', 'L': L_fus, 'R': D_fus/2},
        {'m': m_tail, 'x': X_tail, 'y': 0, 'z': Z_tail, 'type': 'point'},
        {'m': m_gear, 'x': X_gear, 'y': 0, 'z': Z_gear, 'type': 'point'},
        {'m': m_prop, 'x': X_eng, 'y': 0, 'z': Z_eng, 'type': 'sym_point', 'dy': 0.3*b}, # Moteurs à 30% envergure
        {'m': m_sys+m_pay, 'x': X_sys, 'y': 0, 'z': Z_pay, 'type': 'line_x', 'L': 0.8*L_fus} # Réparti
    ]
    
    Ixx_tot, Iyy_tot, Izz_tot, Ixz_tot = Ixx_w, Iyy_w, Izz_w, Ixz_w
    MX_tot, MZ_tot = MX_w, MZ_w
    M_tot = mass_wing_fuel
    
    for comp in components:
        dm = comp['m']
        dx = comp['x']
        dy = comp['y']
        dz = comp['z']
        
        M_tot += dm
        MX_tot += dm * dx
        MZ_tot += dm * dz
        
        # Inerties propres + Transport
        di_xx, di_yy, di_zz = 0, 0, 0
        
        if comp['type'] == 'cyl_x': # Fuselage
            di_xx = dm * (comp['R']**2)
            di_yy = dm * (comp['L']**2 / 12)
            di_zz = dm * (comp['L']**2 / 12)
        elif comp['type'] == 'line_x': # Payload réparti
            di_yy = dm * (comp['L']**2 / 12)
            di_zz = dm * (comp['L']**2 / 12)
        elif comp['type'] == 'sym_point': # Moteurs (+Y et -Y)
            dist_y = comp['dy']
            # Ixx = m(y^2 + z^2) -> 2 * (m/2) * (dist_y^2 + z^2) = m * (dist_y^2 + z^2)
            Ixx_tot += dm * (dist_y**2 + dz**2)
            Iyy_tot += dm * (dx**2 + dz**2)
            Izz_tot += dm * (dx**2 + dist_y**2)
            Ixz_tot += dm * (dx * dz)
            continue # Déjà traité

        # Transport Huygens (Standard)
        Ixx_tot += di_xx + dm * (dy**2 + dz**2)
        Iyy_tot += di_yy + dm * (dx**2 + dz**2)
        Izz_tot += di_zz + dm * (dx**2 + dy**2)
        Ixz_tot += dm * (dx * dz)

    # 4. Calcul CG Global
    X_CG = MX_tot / M_tot
    Z_CG = MZ_tot / M_tot
    
    # 5. Transport au CG Avion (I_CG = I_Nez - M * d^2)
    Ixx_cg = Ixx_tot - M_tot * (Z_CG**2)
    Iyy_cg = Iyy_tot - M_tot * (X_CG**2 + Z_CG**2)
    Izz_cg = Izz_tot - M_tot * (X_CG**2)
    Ixz_cg = Ixz_tot - M_tot * (X_CG * Z_CG)
    
    Ixy_cg = 0.0
    Iyz_cg = 0.0 # Symétrie latérale supposée parfaite
    
    return [Ixx_cg, Iyy_cg, Izz_cg, Ixy_cg, Iyz_cg, Ixz_cg], (X_CG, Z_CG)

# =============================================================================
# MAIN
# =============================================================================

def main():
    os.system('cls' if os.name == 'nt' else 'clear')
    print("="*80)
    print("      INTEGRATION AVION - MISSION LILLE -> COPENHAGUE")
    print("="*80)

    # --- 1. INPUTS ---
    print("\n--- A. CIBLE & MISSION ---")
    mtow_target = input_float(" > MTOW Cible (Structurelle) [kg]", 77000.0)
    n_pax       = input_float(" > Nombre de passagers", 150)
    
    range_km = DIST_LILLE_COPENHAGUE
    print(f" > Trajet : Lille -> Copenhague (~{range_km} km)")
    
    mach_cr     = input_float(" > Mach de croisière", 0.78)
    sfc_eng     = input_float(" > SFC Moteur [mg/N/s]", 16.0)
    lod_cruise  = input_float(" > Finesse de croisière (L/D)", 17.0)

    print("\n--- B. GÉOMÉTRIE AILE ---")
    b_wing      = input_float(" > Envergure [m]", 34.1)
    c_root      = input_float(" > Corde emplanture [m]", 6.0)
    tc_wing     = input_float(" > Épaisseur relative (t/c) [%]", 15.0) / 100.0
    taper       = input_float(" > Effilement (Taper) [0-1]", 0.3)
    sweep       = input_float(" > Flèche (25%) [deg]", 30.0)
    
    # Calcul Surface Ailaire
    S_calc = (b_wing / 2) * c_root * (1 + taper)
    print(f" [INFO] Surface Ailaire Calculée S = {S_calc:.2f} m²")

    print("\n--- C. GÉOMÉTRIE FUSELAGE ---")
    l_fus       = input_float(" > Longueur Fuselage [m]", 37.0)
    d_fus       = input_float(" > Diamètre Fuselage [m]", 4.0)

    mass_payload = n_pax * 95.0

    # --- 2. CALCUL CARBURANT ---
    mass_fuel, trip_fuel = calcul_carburant_breguet(mtow_target, range_km, mach_cr, sfc_eng, lod_cruise)

    # --- 3. OPTIMISATION ---
    print("\n" + "="*80)
    print("   LANCEMENT DES OPTIMISATIONS STRUCTURELLES")
    print("="*80)

    # Aile
    mzfw_est = mtow_target - mass_fuel 
    geo_wing_opt = {'b': b_wing, 'chord_root': c_root, 'tc': tc_wing, 'nz': 3.75}
    print("\n[1/2] Optimisation Aile...")
    m_wing_opt, _, success_wing = wing_opti.optimize_wing_structure(mtow_target, geo_wing_opt)
    
    # Fuselage
    print("\n[2/2] Optimisation Fuselage...")
    m_fus_opt, _, success_fus = fuselage_opti.optimize_fuselage(mtow_target, l_fus, d_fus)

    # Corrections Masses
    m_wing_total = m_wing_opt * 1.35
    m_fus_total  = m_fus_opt * 1.40 

    # --- 4. BILAN DE MASSE DÉTAILLÉ ---
    m_empennage = mtow_target * RATIO_EMPENNAGE
    m_train     = mtow_target * RATIO_TRAIN
    m_nacelles  = mtow_target * RATIO_NACELLES
    m_systemes  = mtow_target * RATIO_SYSTEMES
    m_amenag    = mtow_target * RATIO_AMENAGEMENT
    m_moteurs   = mtow_target * 0.065 

    oew_calc    = m_wing_total + m_fus_total + m_empennage + m_train + m_nacelles + m_moteurs + m_systemes + m_amenag
    mtow_calc   = oew_calc + mass_payload + mass_fuel

    # --- 5. RÉSULTATS GÉOMÉTRIQUES & INERTIELS ---
    m_comps = {
        'wing': m_wing_total, 'fuselage': m_fus_total, 'empennage': m_empennage,
        'train': m_train, 'propulsion': m_moteurs + m_nacelles,
        'systemes': m_systemes + m_amenag, 'payload': mass_payload, 'fuel': mass_fuel
    }
    geo_inertie = {
        'b': b_wing, 'c_root': c_root, 'taper': taper, 'sweep_25': sweep,
        'L_fus': l_fus, 'D_fus': d_fus
    }
    
    I_tensor, CG_pos = calcul_inertie_globale(m_comps, geo_inertie)

    # --- AFFICHAGE FINAL (Tableaux) ---
    os.system('cls' if os.name == 'nt' else 'clear')
    print("="*80)
    print(f"      BILAN FINAL - MISSION LILLE -> COPENHAGUE ({range_km} km)")
    print("="*80)
    
    print("\n1. RÉPARTITION DE MASSE DÉTAILLÉE")
    print("-" * 80)
    print(f"{'COMPOSANT':<25} | {'MASSE (kg)':>12} | {'% MTOW':>10} | {'SOURCE MÉTHODE':<25}")
    print("-" * 80)
    print(f"{'Aile (Structure)':<25} | {m_wing_total:12.1f} | {m_wing_total/mtow_calc*100:9.1f}% | {'Opti. Physique x1.35'}")
    print(f"{'Fuselage (Structure)':<25} | {m_fus_total:12.1f} | {m_fus_total/mtow_calc*100:9.1f}% | {'Opti. Physique x1.40'}")
    print(f"{'Empennages':<25} | {m_empennage:12.1f} | {m_empennage/mtow_calc*100:9.1f}% | {'Statistique'}")
    print(f"{'Train Atterrissage':<25} | {m_train:12.1f} | {m_train/mtow_calc*100:9.1f}% | {'Statistique'}")
    print(f"{'Propulsion (Mot+Nac)':<25} | {(m_moteurs+m_nacelles):12.1f} | {(m_moteurs+m_nacelles)/mtow_calc*100:9.1f}% | {'Statistique'}")
    print(f"{'Systèmes & Aménag.':<25} | {(m_systemes+m_amenag):12.1f} | {(m_systemes+m_amenag)/mtow_calc*100:9.1f}% | {'Statistique'}")
    print("-" * 80)
    print(f"{'OEW (Masse à Vide)':<25} | {oew_calc:12.1f} | {oew_calc/mtow_calc*100:9.1f}% | {'SOUS-TOTAL'}")
    print(f"{'Charge Utile (150 Pax)':<25} | {mass_payload:12.1f} | {mass_payload/mtow_calc*100:9.1f}% | {'Donnée Entrée'}")
    print(f"{'Carburant (Mission+Res)':<25} | {mass_fuel:12.1f} | {mass_fuel/mtow_calc*100:9.1f}% | {'Bréguet (Lille-Copen)'}")
    print("=" * 80)
    print(f"{'MTOW CALCULÉE':<25} | {mtow_calc:12.1f} | {'100%':>10} | {'SOMME TOTALE'}")
    print(f"{'MTOW CIBLE':<25} | {mtow_target:12.1f} | {'':>10} | {'Donnée Entrée'}")
    
    delta = mtow_target - mtow_calc
    print("\n   >>> VERDICT : ", end="")
    if delta >= 0:
        print(f"SUCCÈS ! L'avion est capable d'effectuer la mission (Marge: {delta:.0f} kg)")
    else:
        print(f"ATTENTION ! L'avion est trop lourd pour la structure (Surcharge: {abs(delta):.0f} kg)")

    print("\n2. CARACTÉRISTIQUES POUR SIMULATEUR")
    print("-" * 80)
    print(f"SURFACE AILAIRE (S)   : {S_calc:.2f} m²")
    print(f"CENTRE DE GRAVITÉ (CG): X={CG_pos[0]:.2f}m (depuis le nez), Z={CG_pos[1]:.2f}m")
    print("-" * 40)
    print("TENSEUR D'INERTIE (kg.m²) au CG, Axes Avion (X arrière, Y droite, Z bas)")
    print(f"   Ixx (Roulis)  : {I_tensor[0]:.4E}")
    print(f"   Iyy (Tangage) : {I_tensor[1]:.4E}")
    print(f"   Izz (Lacet)   : {I_tensor[2]:.4E}")
    print(f"   Ixz (Couplage): {I_tensor[5]:.4E}")
    print(f"   Ixy = Iyz     : 0.00 (Symétrie)")
    print("=" * 80)

    # Graphique
    labels = ['Aile', 'Fuselage', 'Autres', 'Propulsion', 'Systèmes', 'Payload', 'Fuel']
    sizes = [m_wing_total, m_fus_total, (m_empennage+m_train), (m_moteurs+m_nacelles), (m_systemes+m_amenag), mass_payload, mass_fuel]
    colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99', '#c2c2f0','#ffb3e6', '#76d7c4']
    
    try:
        plt.figure(figsize=(10, 6))
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
        plt.title(f"Répartition de Masse (Lille-Copenhague) - MTOW={mtow_calc/1000:.1f}t")
        plt.axis('equal')
        plt.show()
    except:
        pass

if __name__ == "__main__":
    main()