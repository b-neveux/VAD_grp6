A = [-0.00255,   6.26,    -0.0486, -9.81;
     -0.000247, -0.488,    0.993,   0;
      0,        -13.3,    -1.50,    0;
      0,         0,        1,       0];

B = [-0.00637,  1.37e-6;   
     -0.00087, -8.46e-10;
     -0.273,    0;
      0,        0];

C = [1, 0, 0, 0;   % Mesure de V (1er état)
     0, 0, 0, 1];  % Mesure de theta (4ème état)
D = [0, 0; 
     0, 0];


disp('Pôles en boucle ouverte :');
disp(eig(A));

Poles_BF = [-0.05 + 0.1j, -0.05 - 0.1j, ...  % Dynamique lente (Vitesse)
            -2.0 + 2.0j,  -2.0 - 2.0j];      % Dynamique rapide (Tangage)

K = place(A, B, Poles_BF);

disp('Gain de retour K (2x4) :');
disp(K);

Poles_Obs = [-0.5, -0.6, -10, -12]; 

L = place(A', C', Poles_Obs)';

disp('Gain de l''observateur L (4x2) :');
disp(L);

A_aug = [A, zeros(4, 2); 
        -C, zeros(2, 2)];

B_aug = [B; zeros(2, 2)];

Poles_Int = [Poles_BF, -0.5, -0.6]; 

K_aug = place(A_aug, B_aug, Poles_Int);

% Séparation des gains
K_p = K_aug(:, 1:4); % Gain proportionnel (sur les états)
K_i = K_aug(:, 5:6); % Gain intégral (sur les erreurs)

disp('Gain Proportionnel K_p :'); disp(K_p);
disp('Gain Intégral K_i :'); disp(K_i);

%Simulation (Réponse à un échelon de Vitesse)
sys_cl = ss(A_aug - B_aug*K_aug, [zeros(4,2); eye(2)], [C, zeros(2,2)], 0);

% Simulation d'une consigne de +10 m/s sur V, 0 sur theta
t = 0:0.1:100;
u_ref = [10; 0]; 
[y, t, x] = lsim(sys_cl, repmat(u_ref', length(t), 1), t);

figure;
subplot(2,1,1); plot(t, y(:,1)); grid on; title('Réponse Vitesse (V)'); ylabel('m/s');
subplot(2,1,2); plot(t, y(:,2)); grid on; title('Réponse Assiette (\theta)'); ylabel('rad');	

