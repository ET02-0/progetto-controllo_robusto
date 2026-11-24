%% ========================================================================
% SETUP H-INFINITY - Elicottero 2DoF
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
Actuators_nom_MIMO = blkdiag(G_act_nom, G_act_nom);
Actuators_unc_MIMO = blkdiag(G_actuator_unc, G_actuator_unc);

% Plant completo
G_nominal = minreal(P_track_nom * Actuators_nom_MIMO, 1e-7);
G_uncertain = P_track_unc * Actuators_unc_MIMO;

disp('Plant nominale e incerto per il tracking creati con successo.');

%% 3. NORMALIZZAZIONE (SCALING)
% Scaliamo il sistema per rendere le variabili adimensionali per l'ottimizzatore
scale_alpha = deg2rad(2);
scale_beta  = deg2rad(3);
scale_F1 = 0.5;
scale_F2 = 0.5;

Sy = diag([scale_alpha, scale_beta]);
Du = diag([scale_F1, scale_F2]);

Dy_inv = inv(Sy);
Du_inv = inv(Du);

% Plant normalizzato
A_reg = A_nom;

P_track_reg = ss(A_reg, B_nom, C_track, D_track);
G_nominal_reg = minreal(P_track_reg * Actuators_nom_MIMO, 1e-7);
G_scaled = minreal(Dy_inv * G_nominal_reg * Du, 1e-7);
G_uncertain_scaled = Dy_inv * G_uncertain * Du;

keepBlocks = {
    'J_alpha'
    'l'
    'omega_n'
    'tau_d'
};
G_unc_mu = G_uncertain;

allBlocks = fieldnames(G_unc_mu.Uncertainty);

for i = 1:numel(allBlocks)
    blk = allBlocks{i};

    if ~ismember(blk, keepBlocks)
        nominalValue = G_unc_mu.Uncertainty.(blk).NominalValue;
        G_unc_mu = usubs(G_unc_mu, blk, nominalValue);
    end
end

G_unc_mu_scaled = Dy_inv * G_unc_mu * Du;

% Pianta con ingressi [u_s; d_s]: d = coppie aerodinamiche (entrano sugli stati, non passano dagli attuatori)
Dd  = diag([aero.alpha.amplitude, aero.beta.amplitude]);       % [N*m]
Pd  = ss(A_nom,[B_nom Bd_nom],C_track,zeros(2,4)) * blkdiag(Actuators_nom_MIMO, eye(2));
Gx  = minreal(Dy_inv*Pd*blkdiag(Du,Dd), 1e-7);                 % [u_s; d_s] -> y_s
Gdd = Gx(:,3:4);                                               % d_s -> y_s
disp(norm(Gx(:,1:2)-G_scaled, inf));                           % controllo: deve essere ~0

disp('Normalizzazione completata.');

%% ------------------------------------------------------------------------
% WS - PESO SULLA SENSIBILITA' S
% -------------------------------------------------------------------------

weight.Ms_alpha = 2.0;
weight.Ms_beta  = 2.0;

weight.As_alpha = 0.001;
weight.As_beta  = 0.001;

weight.wb_alpha = 2.0;
weight.wb_beta  = 1.5;

WS_alpha = ...
    (s/weight.Ms_alpha + weight.wb_alpha) / ...
    (s + weight.wb_alpha*weight.As_alpha);

WS_beta = ...
    (s/weight.Ms_beta + weight.wb_beta) / ...
    (s + weight.wb_beta*weight.As_beta);

WS = blkdiag(WS_alpha, WS_beta);


%% ------------------------------------------------------------------------
% WU - PESO SULLO SFORZO DI CONTROLLO
% -------------------------------------------------------------------------

WU = ss([], [], [], 0.2*eye(2));

%% ------------------------------------------------------------------------
% WT - PESO SULLA SENSIBILITA' COMPLEMENTARE T
% -------------------------------------------------------------------------

weight.Mt_alpha = 1.2;
weight.Mt_beta  = 0.5;

weight.At_alpha = 0.01;
weight.At_beta  = 0.01;

weight.wt_alpha = 10;
weight.wt_beta  = 10;

WT_alpha = ...
    (s + weight.wt_alpha*weight.At_alpha) / ...
    (s/weight.Mt_alpha + weight.wt_alpha);

WT_beta = ...
    (s + weight.wt_beta*weight.At_beta) / ...
    (s/weight.Mt_beta + weight.wt_beta);

WT = blkdiag(WT_alpha, WT_beta);


%% ------------------------------------------------------------------------
% VERIFICA DEI PESI
% -------------------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('PESI COMUNI H-INFINITY / MU-SYNTHESIS\n');
fprintf('============================================================\n');

disp('WS ='); disp(WS);
disp('WU ='); disp(WU);
disp('WT ='); disp(WT);

fprintf('dcgain(WS) =\n');
disp(dcgain(WS));

fprintf('dcgain(WU) =\n');
disp(dcgain(WU));

fprintf('dcgain(WT) =\n');
disp(dcgain(WT));


%% Grafici

omegaWeights = logspace(-2,3,500);

figure('Name','Pesi comuni H-infinity');

subplot(3,1,1);
sigma(inv(WS),omegaWeights);
grid on;
title('W_S^{-1} - limite sulla sensibilita''');

subplot(3,1,2);
sigma(inv(WU),omegaWeights);
grid on;
title('W_U^{-1} - limite sullo sforzo di controllo');

subplot(3,1,3);
sigma(inv(WT),omegaWeights);
grid on;
title('W_T^{-1} - limite sulla sensibilita'' complementare');


%% ------------------------------------------------------------------------
% PLANT GENERALIZZATO COMUNE
% -------------------------------------------------------------------------

% Definiamo i blocchi base con i loro nomi di I/O
G = G_scaled;
G.InputName = {'u1', 'u2'};
G.OutputName = {'y1', 'y2'};

WS.InputName = {'e1', 'e2'};
WS.OutputName = {'zs1', 'zs2'};

WU.InputName = {'u1', 'u2'};
WU.OutputName = {'zu1', 'zu2'};

WT.InputName = {'y1', 'y2'};
WT.OutputName = {'zt1', 'zt2'};

% Definiamo le relazioni di somma e passaggio dei segnali
% r1, r2 sono i riferimenti esterni (esogeni)
% u1, u2 sono i comandi di controllo dal regolatore
blocks = [G; WS; WU; WT];

% Specifichiamo le equazioni di interconnessione
% e = r - y  -->  e1 = r1 - y1,  e2 = r2 - y2
% z_s = WS * e
% z_u = WU * u
% z_t = WT * y
% y_meas = e (l'errore che entra nel controllore K)

% Ingressi del sistema complessivo P_mix: [r1; r2; u1; u2]
% Uscite del sistema complessivo P_mix: [zs1; zs2; zu1; zu2; zt1; zt2; e1; e2]

inputs = {'r1', 'r2', 'u1', 'u2'};
outputs = {'zs1', 'zs2', 'zu1', 'zu2', 'zt1', 'zt2', 'e1', 'e2'};

% Creazione dei blocchi di somma per gli errori
sum1 = sumblk('e1 = r1 - y1');
sum2 = sumblk('e2 = r2 - y2');

% Interconnessione finale con connect
P_mix = connect(G, WS, WU, WT, sum1, sum2, inputs, outputs);

%% 6. VISUALIZZAZIONE E SALVATAGGIO
omegaWeights = logspace(-2,3,500);
figure('Name','H-infinity weighting functions');
subplot(3,1,1); sigma(inv(WS), omegaWeights); grid on; title('W_S^{-1} (Limite sull''errore)');
subplot(3,1,2); sigma(inv(WU), omegaWeights); grid on; title('W_U^{-1} (Limite sugli attuatori)');
subplot(3,1,3); sigma(inv(WT), omegaWeights); grid on; title('W_T^{-1} (Limite per la robustezza)');

save('HINF_workspace.mat', ...
      'G_nominal', ...
      'G_uncertain', ...
      'G_scaled', ...
      'G_uncertain_scaled', ...
      'G_unc_mu', ...
      'G_unc_mu_scaled', ...
      'WS', ...
      'WU', ...
      'WT', ...
      'Sy', ...
      'Du', ...
      'Dy_inv', ...
      'Du_inv', ...
      'omegaWeights', ...
      'Gx', ...
      'Gdd', ...
      'Dd', ...
      'P_mix',...
      'weight');

disp('Setup completato e HINF_workspace.mat salvato!');
