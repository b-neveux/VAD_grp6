clear all; close all; clc;

% 1. Définition du modèle
A = [-0.488,   0.993;
     -13.3,   -1.50];

B = [-8.76e-4;
     -0.273];

C = [1, 0]; 
D = 0;

% Création du système linéaire
sys = ss(A, B, C, D);

poles_BO = eig(A);
disp('Pôles en boucle ouverte:');
disp(poles_BO);

%2. Placement de pôles (Commande u = -Kx)

%Pôles désirés en boucle fermée
%2x plus rapide que la partie réelle
Poles_BF = [-2 + 3j, -2 - 3j]; 

% Calcul du gain K (ou F dans votre notation où F = -K)
K = place(A, B, Poles_BF);

% Calcul du pré-compensateur G (pour avoir un gain statique unitaire, cela assure que la sortie y suit la consigne v
G = inv(C * inv(B*K - A) * B);

disp('Gain de retour K :');
disp(K);
disp('Gain de pré-compensation G :');
disp(G);

% Système en boucle fermée
sys_cl = ss(A - B*K, B*G, C, D);

%3. Observateur de Luenberger

% Les pôles de l'observateur doivent être plus rapides que ceux du système (5 fois par ex)
Poles_Obs = 5 * real(Poles_BF); % Partie réelle seulement pour éviter les oscillations
Poles_Obs = [-10, -11]; % Exemple arbitraire rapide

% Calcul du gain de l'observateur L
L = place(A', C', Poles_Obs)';

disp('Gain de l''observateur L :');
disp(L);

%L'équation de l'observateur est : d(x_est)/dt = A*x_est + B*u + L(y - C*x_est)

%4. Commande Intégrale

% Matrices augmentées
% A_aug = [A, 0; -C, 0]  (taille 3x3)
% B_aug = [B; 0]         (taille 3x1)
A_aug = [A, zeros(2,1); -C, 0];
B_aug = [B; 0];

% Nouveaux pôles désirés (3 pôles car 3 états maintenant)
% On ajoute un pôle pour l'intégrateur
Poles_Int = [Poles_BF, -5]; 

% Calcul du gain global K_aug = [K_p, K_i]
K_aug = place(A_aug, B_aug, Poles_Int);

% Séparation des gains
K_p = K_aug(1:2); % Gain proportionnel
K_i = K_aug(3);   % Gain intégral

disp('Gain augmenté (proportionnel et intégral) :');
disp(K_aug);

% Simulation de la réponse indicielle avec intégrateur
sys_aug = ss(A_aug - B_aug*K_aug, [0;0;1], [C, 0], 0);
figure;
step(sys_aug);
title('Réponse indicielle avec commande intégrale');
grid on;