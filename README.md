# 📘 Projet VAD – Groupe 1

**Analyse et modélisation aérodynamique d’un aéronef**

Ce dépôt contient le travail du Groupe 6 dans le cadre du projet de VAD.
Il centralise le code Python, les modèles aérodynamiques, les données et la documentation technique.

---
## Comment lancer le projet
1. Lancer le programme structure.py dans vad_grp6/src
2. Rentrer les informations utiles à la modélisation
3. Enregistrer et visualiser les résultats dans analysis

## 🔧 Objectifs du projet (rapport intermédiaire)

- Définir les caractéristiques géométriques et massiques d’un aéronef.
- Construire un modèle aérodynamique (coefficients, polaires…).
- Réaliser des analyses via **XFLR5**.
- Mettre en place des scripts Python pour :
  - traiter les données,
  - calculer des coefficients,
  - modéliser le comportement dynamique
  - optimiser la structure d'une aile et d'un fuselage

-> **Les rapports sont disponibles dans VAD_grp6/docs**

---

## 📁 Organisation du dépôt

```text
VAD_grp6/
│
├── src/                         # Scripts Python du projet
│   ├── aircraft_details.py      # Programme de test, non utile dans ce projet
│   ├── fuselage_opti.py         # Programme d'optimisation du fuselage
│   ├── wing_opti.py             # Programme d'optimisation des ailes
│   ├── structure.py             # Programme principal, fait le lien entre wing_opti.py et fuselage_opti.py
│   └── data/
│       ├── Data_k_bh.py         # Fichier de coefficients k,bh
│       └── Aircraft_Data.xlsx   # Tableau des caractéristiques de plusieurs avions
│   └── analysis/
│       └── plane_name.txt       # Informations issues de l'analyse d'un avion
│
├── models/
│   └── XFLR5_tests.xfl          # Modèles XFLR5
│
├── matlab/
│   ├── VAD_modele_1.m           # Fichier matlab correspondant au modèle 1
│   └── VAD_modele_2.m           # Fichier matlab correspondant au modèle 2
│
├── plane_score/
│   └── flight_dynamics_launcher_2.exe
│
├── docs/
│   ├── VAD_rapport_intermediaire.pdf
│   ├── VAD_rapport_final.pdf
│   └── cours/                   # Cours de VAD
│
└── README.md                    # Documentation du projet