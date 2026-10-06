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

%% 4. DEFINIZIONE DEI PESI FREQUENZIALI (W_S, W_U, W_T)
%% ========================================================================
%  PESI DI PRESTAZIONE PER LA SINTESI MIXED-SENSITIVITY
%
%  Problema:   min || [ W1*S ; W2*K*S ; W3*T ] ||_inf  < 1
%
%  W1 -> forma S: banda, errore a regime, picco di sensitivita (= margini)
%  W2 -> forma K*S: sforzo di controllo e amplificazione del rumore
%  W3 -> forma T: roll-off, robustezza a dinamica non modellata
% ========================================================================

s = tf('s');

%% ------------------------------------------------------------------ WS %%
%  WS_i = (s/M + wb_i)/(s + wb_i*A).
M_S = [1.8; 1.7];
A_S = 0.10;
wb  = [2.0;1.5];

WS_a = (s/M_S(1) + wb(1))/(s + wb(1)*A_S);
WS_b = ((s/M_S(2) + wb(2))/(s + wb(2)*A_S))^2;
WS   = blkdiag(WS_a, WS_b);



%% ------------------------------------------------------------------ WU %%
%  Vincolo: |K*S| <= u_max/e_max in bassa frequenza, con forte penalizzazione
%  in alta frequenza (roll-off del controllore -> rumore IMU non amplificato).

kU = 0.1;
wu = 8;                 % [rad/s] inizio della penalizzazione
Mu = 7;                 % LF: 0.1, HF: 0.7
WU_ch = kU*(s/wu + 1)/(s/(Mu*wu) + 1);
WU = blkdiag(WU_ch, WU_ch);

%% ------------------------------------------------------------------ WT %%
%  WT_i = (s + wbt/M_T)/(A_T*s + wbt)
M_T = 1.50;
A_T = 0.01;
wbt = [22; 18];                     % [rad/s]

WT = blkdiag( (s + wbt(1)/M_T)/(A_T*s + wbt(1)), ...
              (s + wbt(2)/M_T)/(A_T*s + wbt(2)) );
WT.u = {'y(1)','y(2)'};  
WT.y = {'z3(1)','z3(2)'};

%% -------------------------------------------------------- visualizzazione
w = logspace(-3,3,800);
figure('Name','Pesi di prestazione');
sigma(1/WS(1,1),'b', 1/WU(1,1),'r', 1/WT(1,1),'g', w); grid on
legend('1/W1 (vincolo su S)','1/W2 (vincolo su KS)','1/W3 (vincolo su T)', ...
       'Location','best');
title('Vincoli di progetto (canale pitch)');

fprintf(['Specifiche: wb = [%.1f %.1f] rad/s, ||S||inf <= %.1f, ' ...
         '||T||inf <= %.1f, roll-off T = %.0f dB\n'], wb(1),wb(2),M_S,M_T, ...
         20*log10(A_T));
fprintf('=== s02 completato: ora esegui s03_mixsyn ===\n\n');

%% ------------------------------------------------------------ COME TARARE
%  Se gamma > 1 in s03/s04:
%    1. abbassa wb (prima cosa da provare: wb = [2 1.5]);
%    2. alza M_S a 2.5 e M_T a 2.5;
%    3. alza A_T (roll-off meno aggressivo, es. 0.05);
%    4. alza kU (piu autorita al controllore) SOLO se le saturazioni lo
%       consentono: verificalo poi in s07 sul modello non lineare.
%  Se gamma << 1 le specifiche sono poco ambiziose: alza wb.

%% 5. COSTRUZIONE DEL PLANT GENERALIZZATO (P_mix)
P_mix = augw(G_scaled, WS, WU, WT);             

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
      'omegaWeights');

disp('Setup completato e HINF_workspace.mat salvato!');
