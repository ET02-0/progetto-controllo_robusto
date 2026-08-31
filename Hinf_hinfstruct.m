%{
%% ========================================================================
% SINTESI H-INFINITY STRUTTURATA (HINFSTRUCT)
% ========================================================================
disp('============================================================');
disp(' SINTESI HINFSTRUCT - PID + LEAD/LAG MIMO');
disp('============================================================');

load('HINF_workspace.mat');
s = tf('s');
I2 = eye(2);

%% 1. INIZIALIZZAZIONE PID (TUNABLE)
% Fattori di conversione per inizializzare i PID con valori ragionevoli
alphaPhysToScaled = scale_alpha / scale_F1;
betaPhysToScaled  = scale_beta / scale_F2;

Kalpha_tunable = tunablePID('Kalpha', 'PID');
Kalpha_tunable.Kp.Value = alphaPhysToScaled * (-0.2); 
Kalpha_tunable.Ki.Value = alphaPhysToScaled * (1.4);
Kalpha_tunable.Tf.Value = 0.02;
Kalpha_tunable.Ki.Minimum = 1e-10; % Evita integratori nulli
Kalpha_tunable.Tf.Minimum = 1e-4;

Kbeta_tunable = tunablePID('Kbeta', 'PID');
Kbeta_tunable.Kp.Value = betaPhysToScaled * (0.2);
Kbeta_tunable.Ki.Value = betaPhysToScaled * (0.04);
Kbeta_tunable.Tf.Value = 0.02;
Kbeta_tunable.Ki.Minimum = 1e-10;
Kbeta_tunable.Tf.Minimum = 1e-4;

%% 2. COMPENSATORI DINAMICI (LEAD/LAG)
Falpha0 = (1 + s/2.5) / (1 + s/8);
Fbeta0  = (1 + s/2.0) / (1 + s/6);

Falpha_tunable = tunableTF('Falpha', tf(Falpha0));
Fbeta_tunable  = tunableTF('Fbeta',  tf(Fbeta0));

% Assemblaggio dei blocchi diagonali (Canali indipendenti)
Calpha_tunable = Falpha_tunable * Kalpha_tunable;
Cbeta_tunable  = Fbeta_tunable  * Kbeta_tunable;
Kdiag_tunable  = blkdiag(Calpha_tunable, Cbeta_tunable);

%% 3. DISACCOPPIATORE STATICO (DECOUPLER)
Gdc_scaled = dcgain(G_scaled);
Ddec_tunable = tunableGain('Ddec', 2, 2);
Ddec_tunable.Gain.Value = inv(Gdc_scaled);

% Controllore strutturato completo (Decoupler + PID_LeadLag)
K_pidcomp_tunable = Ddec_tunable * Kdiag_tunable;

%% 4. COSTRUZIONE DEL LOOP PARAMETRICO DA OTTIMIZZARE
L_tunable  = G_scaled * K_pidcomp_tunable;
S_tunable  = feedback(I2, L_tunable);
T_tunable  = feedback(L_tunable, I2);
KS_tunable = K_pidcomp_tunable * S_tunable;

% Vettore pesato per hinfstruct
CL_tunable = [WS * S_tunable; 
              WU * KS_tunable; 
              WT * T_tunable];

%% 5. OTTIMIZZAZIONE
% Usiamo RandomStart per evitare minimi locali
optsHinfStruct = hinfstructOptions('Display', 'final', 'RandomStart', 15);
[CL_tuned, gamma_struct, info_struct] = hinfstruct(CL_tunable, optsHinfStruct);

%% 6. ESTRAZIONE E RICONVERSIONE
% Estraiamo i blocchi ottimizzati
Ddec_tuned   = getBlockValue(CL_tuned, 'Ddec');
Kalpha_tuned = getBlockValue(CL_tuned, 'Kalpha');
Kbeta_tuned  = getBlockValue(CL_tuned, 'Kbeta');
Falpha_tuned = getBlockValue(CL_tuned, 'Falpha');
Fbeta_tuned  = getBlockValue(CL_tuned, 'Fbeta');

% Creiamo il controllore lineare finale scalato
Calpha_scaled = minreal(ss(Falpha_tuned) * ss(Kalpha_tuned), 1e-7);
Cbeta_scaled  = minreal(ss(Fbeta_tuned)  * ss(Kbeta_tuned),  1e-7);
K_struct_scaled = minreal(ss(Ddec_tuned) * blkdiag(Calpha_scaled, Cbeta_scaled), 1e-7);

% Riportiamo alle unità fisiche
K_struct = minreal(Du * K_struct_scaled * Dy_inv, 1e-7);

Ddec_phys = inv(dcgain(G_scaled));

fprintf('\nGamma hinfstruct = %.4f\n', gamma_struct);
fprintf('Ordine del controllore K_struct = %d\n', order(K_struct_scaled));

save('HINF_controllers_struct.mat', 'K_struct', 'K_struct_scaled', 'gamma_struct');
disp('Sintesi Strutturata completata!');
%}
%% ========================================================================
% SINTESI HINFSTRUCT CORRETTA - NO IMPROPER PLANT ERROR
% ========================================================================
disp('============================================================');
disp(' SINTESI HINFSTRUCT - PID + LEAD/LAG MIMO');
disp('============================================================');

load('HINF_workspace.mat');
s = tf('s');
I2 = eye(2);

%% 1. PID PARAMETRICI (CONVERTITI IN SS PARAMETRICO)
Kalpha_pid = tunablePID('Kalpha', 'PID');
Kalpha_pid.Kp.Value = 1.0;  Kalpha_pid.Kp.Minimum = 0;
Kalpha_pid.Ki.Value = 0.5;  Kalpha_pid.Ki.Minimum = 1e-4;
Kalpha_pid.Kd.Value = 0.1;  Kalpha_pid.Kd.Minimum = 0;
Kalpha_pid.Tf.Value = 0.02; Kalpha_pid.Tf.Minimum = 1e-3;

Kbeta_pid = tunablePID('Kbeta', 'PID');
Kbeta_pid.Kp.Value = 2.0;   Kbeta_pid.Kp.Minimum = 0;
Kbeta_pid.Ki.Value = 0.5;   Kbeta_pid.Ki.Minimum = 1e-4;
Kbeta_pid.Kd.Value = 0.2;   Kbeta_pid.Kd.Minimum = 0;
Kbeta_pid.Tf.Value = 0.02;  Kbeta_pid.Tf.Minimum = 1e-3;

%% 2. LEAD/LAG CON TUNABLESS (ORDINE 1 INIZIALIZZATO)
Falpha_init = (1 + s/3.0) / (1 + s/12.0);
Fbeta_init  = (1 + s/2.0) / (1 + s/10.0);

% Inizializzazione corretta del blocco tunableSS
Falpha_tunable = tunableSS('Falpha', ss(Falpha_init));
Fbeta_tunable  = tunableSS('Fbeta',  ss(Fbeta_init));

% Cascata canali indipendenti
Calpha_tunable = Falpha_tunable * Kalpha_pid;
Cbeta_tunable  = Fbeta_tunable  * Kbeta_pid;
Kdiag_tunable  = blkdiag(Calpha_tunable, Cbeta_tunable);

%% 3. DISACCOPPIATORE STATICO
D_phys_init = inv(B_nom([2, 4], :));
D_scaled_init = Du_inv * D_phys_init * Dy;
Ddec_tunable = tunableGain('Ddec', D_scaled_init);

K_pidcomp_tunable = Ddec_tunable * Kdiag_tunable;

%% 4. COSTRUZIONE INTERCONNESSIONE CLOSED-LOOP
L_tunable  = G_scaled * K_pidcomp_tunable;
S_tunable  = feedback(I2, L_tunable);
T_tunable  = feedback(L_tunable, I2);
KS_tunable = K_pidcomp_tunable * S_tunable;

CL_tunable = [WS * S_tunable; 
              WU * KS_tunable; 
              WT * T_tunable];

%% 5. OTTIMIZZAZIONE HINFSTRUCT
opts = hinfstructOptions('Display', 'final', 'RandomStart', 10, 'MaxIter', 300);
[CL_tuned, gamma_struct] = hinfstruct(CL_tunable, opts);

%% 6. ESTRAZIONE VALORI PER SIMULINK
Ddec_tuned   = getBlockValue(CL_tuned, 'Ddec');
Kalpha_tuned = getBlockValue(CL_tuned, 'Kalpha');
Kbeta_tuned  = getBlockValue(CL_tuned, 'Kbeta');
Falpha_ss    = getBlockValue(CL_tuned, 'Falpha');
Fbeta_ss     = getBlockValue(CL_tuned, 'Fbeta');

% Conversione Transfer Fcn per i blocchi di Simulink
Falpha_tf = tf(Falpha_ss);
Fbeta_tf  = tf(Fbeta_ss);

% Matrice disaccoppiatore per Simulink
Ddec_phys = Du * Ddec_tuned * Dy_inv;

fprintf('\n=== SINTESI COMPLETATA CON SUCCESSO ===\n');
fprintf('Gamma H-inf ottenuto: %.4f\n', gamma_struct);