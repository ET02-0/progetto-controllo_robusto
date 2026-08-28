%% ========================================================================
% ANALISI NOMINALE H-INFINITY E CONFRONTO CONTROLLORI
% ========================================================================
disp('============================================================');
disp(' ANALISI PRESTAZIONI: TEMPO E FREQUENZA');
disp('============================================================');

% 1. Caricamento Dati
load('HINF_workspace.mat');
load('HINF_controllers_full.mat');
load('HINF_controllers_struct.mat');

I2 = eye(2);
tStep = 0:0.01:15; % Vettore di tempo per la simulazione

% Organizzazione dei controllori
controllersScaled = {K_mix_scaled, K_hinfsyn_scaled, K_struct_scaled};
controllersPhysical = {K_mix, K_hinfsyn, K_struct};
controllerNames = {'Mixsyn (Full)', 'Hinfsyn (Full)', 'Hinfstruct (PID)'};
nContr = length(controllerNames);

% Array per salvare le metriche temporali
pitchSettling = zeros(nContr,1); pitchOvershoot = zeros(nContr,1);
yawSettling = zeros(nContr,1);   yawOvershoot = zeros(nContr,1);

% Figure per le risposte al gradino
fig_pitch = figure('Name','H-infinity - Pitch Tracking'); hold on; grid on;
fig_yaw   = figure('Name','H-infinity - Yaw Tracking'); hold on; grid on;

% Figure per l'analisi in frequenza
fig_S  = figure('Name','Funzione di Sensibilità S'); hold on; grid on;
fig_KS = figure('Name','Sensibilità del Controllo KS'); hold on; grid on;
fig_T  = figure('Name','Sensibilità Complementare T'); hold on; grid on;

%% 2. ITERAZIONE SUI CONTROLLORI E CALCOLO
for k = 1:nContr
    Ks = controllersScaled{k};
    Kp = controllersPhysical{k};
    
    % --- Analisi in Frequenza (Modello Scalato) ---
    Ls = G_scaled * Ks;
    Ss = minreal(feedback(I2, Ls), 1e-7);
    Ts = minreal(feedback(Ls, I2), 1e-7);
    KSs = minreal(Ks * Ss, 1e-7);
    
    % Plot Valori Singolari
    figure(fig_S);  sigma(Ss); 
    figure(fig_KS); sigma(KSs);
    figure(fig_T);  sigma(Ts);
    
    % --- Analisi Temporale (Modello Fisico) ---
    Lp = G_nominal * Kp;
    Tp = minreal(feedback(Lp, I2), 1e-7);
    
    % Gradino su Pitch (alpha)
    [y_pitch, ~] = step(Tp(1,1) * scale_alpha, tStep);
    figure(fig_pitch); plot(tStep, rad2deg(y_pitch), 'LineWidth', 1.5);
    info_p = stepinfo(y_pitch, tStep, scale_alpha, 'SettlingTimeThreshold', 0.02);
    pitchSettling(k) = info_p.SettlingTime;
    pitchOvershoot(k) = info_p.Overshoot;
    
    % Gradino su Yaw (beta)
    [y_yaw, ~] = step(Tp(2,2) * scale_beta, tStep);
    figure(fig_yaw); plot(tStep, rad2deg(y_yaw), 'LineWidth', 1.5);
    info_y = stepinfo(y_yaw, tStep, scale_beta, 'SettlingTimeThreshold', 0.02);
    yawSettling(k) = info_y.SettlingTime;
    yawOvershoot(k) = info_y.Overshoot;
end

%% 3. COMPLETAMENTO GRAFICI FREQUENZIALI
omegaHinf = logspace(-2, 3, 500);

figure(fig_S); 
sigma(inv(WS), omegaHinf, 'k--'); 
title('Sensibilità S(j\omega)'); legend([controllerNames, 'W_S^{-1}']);

figure(fig_KS); 
sigma(inv(WU), omegaHinf, 'k--'); 
title('Sensibilità Controllo KS(j\omega)'); legend([controllerNames, 'W_U^{-1}']);

figure(fig_T); 
sigma(inv(WT), omegaHinf, 'k--'); 
title('Sensibilità Complementare T(j\omega)'); legend([controllerNames, 'W_T^{-1}']);

%% 4. COMPLETAMENTO GRAFICI TEMPORALI E STAMPA METRICHE
figure(fig_pitch);
yline(rad2deg(scale_alpha), 'k--', 'Riferimento');
title('Inseguimento Gradino Pitch'); xlabel('Tempo [s]'); ylabel('Pitch [deg]');
legend(controllerNames, 'Location', 'best');

figure(fig_yaw);
yline(rad2deg(scale_beta), 'k--', 'Riferimento');
title('Inseguimento Gradino Yaw'); xlabel('Tempo [s]'); ylabel('Yaw [deg]');
legend(controllerNames, 'Location', 'best');

disp('--- METRICHE TEMPORALI ---');
for k = 1:nContr
    fprintf('\n%s:\n', controllerNames{k});
    fprintf('  Pitch -> Ts = %.2f s, Overshoot = %.2f %%\n', pitchSettling(k), pitchOvershoot(k));
    fprintf('  Yaw   -> Ts = %.2f s, Overshoot = %.2f %%\n', yawSettling(k), yawOvershoot(k));
end