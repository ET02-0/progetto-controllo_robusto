%% ========================================================================
% SINTESI H-INFINITY FULL-ORDER
% ========================================================================
disp('============================================================');
disp(' SINTESI H-INFINITY: MIXSYN E HINFSYN');
disp('============================================================');

% Carichiamo il setup generato in precedenza
load('HINF_workspace.mat');

% Opzioni per visualizzare i dettagli dell'ottimizzazione
optsHinf = hinfsynOptions('Display','on');

%% 1. SINTESI CON MIXSYN
disp('--- Esecuzione MIXSYN ---');
% mixsyn prende direttamente il plant scalato e i tre pesi
[K_mix_scaled, CL_mix, gamma_mix, info_mix] = mixsyn(G_scaled, WS, WU, WT, optsHinf);

% Pulizia numerica degli stati non necessari
K_mix_scaled = minreal(ss(K_mix_scaled), 1e-7);

% Riportiamo il controllore alle unità di misura fisiche dell'elicottero
K_mix = minreal(Du * K_mix_scaled * Dy_inv, 1e-7);

fprintf('Gamma mixsyn = %.4f\n', gamma_mix);
fprintf('Ordine del controllore K_mix = %d\n\n', order(K_mix_scaled));

%% 2. SINTESI CON HINFSYN
disp('--- Esecuzione HINFSYN ---');
% hinfsyn richiede il plant generalizzato P_mix già assemblato
nmeas = 2; % Numero di misure (alpha e beta)
ncont = 2; % Numero di comandi (F1 e F2)

[K_hinfsyn_scaled, CL_hinfsyn, gamma_hinfsyn, info_hinfsyn] = hinfsyn(P_mix, nmeas, ncont, optsHinf);

% Pulizia numerica
K_hinfsyn_scaled = minreal(ss(K_hinfsyn_scaled), 1e-7);

% Riconversione fisica
K_hinfsyn = minreal(Du * K_hinfsyn_scaled * Dy_inv, 1e-7);

fprintf('Gamma hinfsyn = %.4f\n', gamma_hinfsyn);
fprintf('Ordine del controllore K_hinfsyn = %d\n', order(K_hinfsyn_scaled));

%% 3. SALVATAGGIO
save('HINF_controllers_full.mat', 'K_mix', 'K_mix_scaled', 'gamma_mix', ...
     'K_hinfsyn', 'K_hinfsyn_scaled', 'gamma_hinfsyn');
disp('Sintesi completata e controllori salvati!');