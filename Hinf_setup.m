%% ========================================================================
% SETUP H-INFINITY - Elicottero 2DoF
% Adattato per le variabili del workspace di dataset_elicottero
% ========================================================================
disp('============================================================');
disp(' INIZIALIZZAZIONE SETUP H-INFINITY');
disp('============================================================');

%% 1. DEFINIZIONE DEL PLANT PER IL TRACKING (Estrazione alpha e beta)
% Creiamo le matrici di uscita per estrarre solo gli angoli
C_track = [
    1 0 0 0; % Estrae alpha
    0 0 1 0  % Estrae beta
];
D_track = zeros(2,2);

% Plant nominale per il tracking
P_track_nom = ss(A_nom, B_nom, C_track, D_track);
P_track_nom.InputName = {'F1', 'F2'};
P_track_nom.OutputName = {'alpha', 'beta'};

% Plant incerto per il tracking
P_track_unc = uss(A_unc, B_unc, C_track, D_track);
P_track_unc.InputName = {'F1', 'F2'};
P_track_unc.OutputName = {'alpha', 'beta'};

%% 2. INCLUSIONE DINAMICA ATTUATORI
% Usiamo le tue funzioni di trasferimento G_act_nom e G_actuator_unc
% e le assembliamo in una matrice MIMO diagonale (2 ingressi, 2 uscite)
Actuators_nom_MIMO = blkdiag(G_act_nom, G_act_nom);
Actuators_unc_MIMO = blkdiag(G_actuator_unc, G_actuator_unc);

% Plant completo: Attuatori in serie alla dinamica dell'elicottero
G_nominal = minreal(P_track_nom * Actuators_nom_MIMO, 1e-7);
G_uncertain = P_track_unc * Actuators_unc_MIMO;

disp('Plant nominale e incerto per il tracking creati con successo.');

%% 3. NORMALIZZAZIONE (SCALING)
% Scaliamo il sistema per rendere le variabili adimensionali per l'ottimizzatore
scale_alpha = deg2rad(2);
scale_beta  = deg2rad(3);
scale_F1 = 2.5;
scale_F2 = 2.5;

Dy = diag([scale_alpha, scale_beta]);
Du = diag([scale_F1, scale_F2]);

Dy_inv = inv(Dy);
Du_inv = inv(Du);

% Plant normalizzato (Quello su cui H-inf farà i calcoli)
% Creiamo un plant regolarizzato per eliminare il polo quasi-nullo
A_reg = A_nom;
% Aggiungiamo smorzamento sull'asse di yaw (velocità di imbardata, stato 4)
A_reg(4, 4) = A_reg(4, 4) - 0.5; 

P_track_reg = ss(A_reg, B_nom, C_track, D_track);
G_nominal_reg = minreal(P_track_reg * Actuators_nom_MIMO, 1e-7);
G_scaled = minreal(Dy_inv * G_nominal_reg * Du, 1e-7);
G_uncertain_scaled = Dy_inv * G_uncertain * Du;

disp('Normalizzazione completata.');

%% 4. DEFINIZIONE DEI PESI FREQUENZIALI (W_S, W_U, W_T)
s = tf('s');

% --- W_S: Peso sulla Funzione di Sensibilità (Errore di tracking) ---
Ms_alpha = 1.45;  wb_alpha = 3.7;  As_alpha = 0.015;
Ms_beta  = 1.50;  wb_beta  = 5;  As_beta  = 0.015;

WS_alpha = (s/Ms_alpha + wb_alpha) / (s + wb_alpha*As_alpha);
WS_beta  = (s/Ms_beta  + wb_beta)  / (s + wb_beta*As_beta);
WS = blkdiag(WS_alpha, WS_beta);

% --- W_U: Peso sullo Sforzo di Controllo ---
% Lasciamo un margine alto (0.90) per permettere rotori reattivi
WU = ss([], [], [], 0.80*eye(2));

% --- W_T: Peso sulla Sensibilità Complementare (Robustezza al rumore/ritardo) ---
Mt_alpha = 1.50;  wt_alpha = 22;  At_alpha = 0.01;
Mt_beta  = 1.50;  wt_beta  = 18;  At_beta  = 0.01;

WT_alpha = (s + wt_alpha*At_alpha) / (s/Mt_alpha + wt_alpha);
WT_beta  = (s + wt_beta*At_beta)   / (s/Mt_beta  + wt_beta);
WT = blkdiag(WT_alpha, WT_beta);

%% 5. COSTRUZIONE DEL PLANT GENERALIZZATO (P_mix)
% Unisce il plant normalizzato e i pesi nel formato standard H-infinity
P_mix = augw(G_scaled, WS, WU, WT);

%% 6. VISUALIZZAZIONE E SALVATAGGIO
omegaWeights = logspace(-2,3,500);
figure('Name','H-infinity weighting functions');
subplot(3,1,1); sigma(inv(WS), omegaWeights); grid on; title('W_S^{-1} (Limite sull''errore)');
subplot(3,1,2); sigma(inv(WU), omegaWeights); grid on; title('W_U^{-1} (Limite sugli attuatori)');
subplot(3,1,3); sigma(inv(WT), omegaWeights); grid on; title('W_T^{-1} (Limite per la robustezza)');

% Salviamo il workspace per gli script di sintesi e analisi
save('HINF_workspace.mat', 'G_nominal', 'G_uncertain', 'G_scaled', 'G_uncertain_scaled', ...
     'P_mix', 'WS', 'WU', 'WT', 'Dy', 'Du', 'Dy_inv', 'Du_inv');

disp('Setup completato e HINF_workspace.mat salvato!');