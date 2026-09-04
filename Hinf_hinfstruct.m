%% ========================================================================
% H-INFINITY STRUTTURATO - HINFSTRUCT
% ELICOTTERO 2DoF
%
% Struttura:
%   PID_alpha + Lead/Lag_alpha
%   PID_beta  + Lead/Lag_beta
%   MIMO static decoupler
%
% Sintesi mixed-sensitivity:
%
%       [ WS*S  ]
% CL =  [ WU*KS ]
%       [ WT*T  ]
%
% Il progetto viene eseguito sul plant SCALATO G_scaled.
%
% Conversione finale:
%
%       K_phys = Du * K_scaled * Dy_inv
%
% ========================================================================
clc;

disp(' ');
disp('============================================================');
disp('       SINTESI H-INFINITY STRUTTURATA - HINFSTRUCT');
disp('       ELICOTTERO 2DoF');
disp('============================================================');


%% ========================================================================
% 0. CARICAMENTO WORKSPACE
% ========================================================================

load('HINF_workspace.mat');

s  = tf('s');
I2 = eye(2);

fprintf('\nWorkspace caricato.\n');

%% Aggiunta test
%% DIAGNOSTICA PLANT E PESI

fprintf('\n============================================================\n');
fprintf(' DIAGNOSTICA INTEGRATORI PLANT E PESI\n');
fprintf('============================================================\n');

disp('--- G_scaled ---');
disp(zpk(G_scaled));

disp('--- WS ---');
disp(zpk(WS));

disp('--- WU ---');
disp(zpk(WU));

disp('--- WT ---');
disp(zpk(WT));

fprintf('\nStabilità:\n');
fprintf('G_scaled = %d\n',isstable(G_scaled));
fprintf('WS       = %d\n',isstable(WS));
fprintf('WU       = %d\n',isstable(WU));
fprintf('WT       = %d\n',isstable(WT));

pG  = pole(G_scaled);
pWS = pole(WS);
pWU = pole(WU);
pWT = pole(WT);

fprintf('\nNumero poli esattamente/vicino a zero:\n');
fprintf('G_scaled = %d\n',sum(abs(pG)  < 1e-8));
fprintf('WS       = %d\n',sum(abs(pWS) < 1e-8));
fprintf('WU       = %d\n',sum(abs(pWU) < 1e-8));
fprintf('WT       = %d\n',sum(abs(pWT) < 1e-8));
%% Fine aggiunta

%% ========================================================================
% 1. CONTROLLI PRELIMINARI SUL PLANT
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 1. CONTROLLO PLANT SCALATO\n');
fprintf('============================================================\n');

fprintf('Dimensione G_scaled: %d x %d\n', ...
    size(G_scaled,1), size(G_scaled,2));

if ~isequal(size(G_scaled),[2 2])
    error('G_scaled deve essere un plant 2x2.');
end

fprintf('\nDC gain G_scaled:\n');
disp(dcgain(G_scaled));

Gdc_scaled = dcgain(G_scaled);

if any(~isfinite(Gdc_scaled(:)))
    error('dcgain(G_scaled) contiene NaN o Inf.');
end

if rcond(Gdc_scaled) < 1e-10
    error(['dcgain(G_scaled) è numericamente quasi singolare. ' ...
           'Il decoupler statico non può essere inizializzato in modo affidabile.']);
end

fprintf('cond(Gdc_scaled) = %.3e\n',cond(Gdc_scaled));
fprintf('rcond(Gdc_scaled) = %.3e\n',rcond(Gdc_scaled));


%% ========================================================================
% 2. PID ALPHA
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 2. PID ALPHA\n');
fprintf('============================================================\n');

Kalpha_tunable = tunablePID('Kalpha','PID');

% ------------------------------------------------------------------------
% Valori iniziali
%
% ATTENZIONE:
% questi sono valori nel MODELLO SCALATO, perché il controllore è
% inserito direttamente nel loop con G_scaled.
% ------------------------------------------------------------------------

Kalpha_tunable.Kp.Value = -0.215;
Kalpha_tunable.Ki.Value =  0;
Kalpha_tunable.Kd.Value =  0.000;
Kalpha_tunable.Tf.Value =  0.020;

% Vincoli
%Kalpha_tunable.Ki.Minimum = 1e-10;
Kalpha_tunable.Kd.Minimum = 0;
Kalpha_tunable.Tf.Minimum = 1e-4;


%% ========================================================================
% 3. PID BETA
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 3. PID BETA\n');
fprintf('============================================================\n');

Kbeta_tunable = tunablePID('Kbeta','PID');

Kbeta_tunable.Kp.Value = 0.215;
Kbeta_tunable.Ki.Value = 0;
Kbeta_tunable.Kd.Value = 0.000;
Kbeta_tunable.Tf.Value = 0.020;

% Vincoli
%Kbeta_tunable.Ki.Minimum = 1e-10;
Kbeta_tunable.Kd.Minimum = 0;
Kbeta_tunable.Tf.Minimum = 1e-4;


%% ========================================================================
% 4. LEAD/LAG ALPHA
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 4. LEAD/LAG ALPHA\n');
fprintf('============================================================\n');

%       1 + s/z
% F = -----------
%       1 + s/p

Falpha0 = (1 + s/3.0) / (1 + s/12.0);

Falpha_tunable = tunableTF( ...
    'Falpha', ...
    tf(Falpha0));


%% ========================================================================
% 5. LEAD/LAG BETA
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 5. LEAD/LAG BETA\n');
fprintf('============================================================\n');

Fbeta0 = (1 + s/2.0) / (1 + s/10.0);

Fbeta_tunable = tunableTF( ...
    'Fbeta', ...
    tf(Fbeta0));


%% ========================================================================
% 6. COSTRUZIONE DEL CONTROLLORE DIAGONALE
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 6. CONTROLLORE DIAGONALE\n');
fprintf('============================================================\n');

Calpha_tunable = Falpha_tunable * Kalpha_tunable;
Cbeta_tunable  = Fbeta_tunable  * Kbeta_tunable;

Kdiag_tunable = blkdiag( ...
    Calpha_tunable, ...
    Cbeta_tunable);


%% ========================================================================
% 7. DECOUPLER STATICO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 7. DECOUPLER STATICO\n');
fprintf('============================================================\n');

% Inizializzazione:
%
%       Ddec = Gdc_scaled^(-1)
%
% ma calcolata numericamente tramite "\"
%
%       Gdc_scaled * Ddec = I

Ddec_init = Gdc_scaled \ I2;

fprintf('Ddec iniziale (scaled):\n');
disp(Ddec_init);

fprintf('Controllo Gdc_scaled * Ddec_init:\n');
disp(Gdc_scaled * Ddec_init);

Ddec_tunable = tunableGain('Ddec',Ddec_init);


%% ========================================================================
% 8. CONTROLLORE MIMO COMPLETO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 8. CONTROLLORE MIMO COMPLETO\n');
fprintf('============================================================\n');

K_pidcomp_tunable = ...
    Ddec_tunable * Kdiag_tunable;


%% ========================================================================
% 9. LOOP OPEN E CLOSED LOOP PARAMETRICO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 9. COSTRUZIONE CLOSED LOOP\n');
fprintf('============================================================\n');

L_tunable = ...
    G_scaled * K_pidcomp_tunable;

S_tunable = ...
    feedback(I2,L_tunable);

T_tunable = ...
    feedback(L_tunable,I2);

KS_tunable = ...
    K_pidcomp_tunable * S_tunable;


%% ========================================================================
% 10. MIXED SENSITIVITY
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 10. MIXED SENSITIVITY\n');
fprintf('============================================================\n');

% Obiettivo H-infinity:
%
%       WS*S
%       WU*K*S
%       WT*T

CL_tunable = [ ...
    WS*S_tunable;
    WU*KS_tunable;
    WT*T_tunable];

fprintf('\n============================================================\n');
fprintf(' BLOCCHI TUNABLE\n');
fprintf('============================================================\n');

showTunable(CL_tunable);

%% ========================================================================
% 11. OPZIONI HINFSTRUCT
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 11. AVVIO HINFSTRUCT\n');
fprintf('============================================================\n');

opts = hinfstructOptions( ...
    'Display','final', ...
    'RandomStart',15, ...
    'MaxIter',300, ...
    'MinDecay', 1e-5);


%% ========================================================================
% 12. SINTESI H-INFINITY
% ========================================================================

[CL_tuned,gamma_struct,info_struct] = ...
    hinfstruct(CL_tunable,opts);


%% ========================================================================
% 13. RISULTATO DELLA SINTESI
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 13. RISULTATO HINFSTRUCT\n');
fprintf('============================================================\n');

fprintf('Gamma H-inf = %.10f\n',gamma_struct);


%% ========================================================================
% 14. VERIFICA DIRETTA DEL MODELLO RESTITUITO DA HINFSTRUCT
%
% QUESTO BLOCCO È MOLTO IMPORTANTE.
%
% Prima di ricostruire manualmente K, verifichiamo direttamente CL_tuned.
% In questo modo distinguiamo:
%
%   A) problema nella sintesi HINFSTRUCT
%
% da
%
%   B) problema nella successiva estrazione/ricostruzione del controllore.
%
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 14. VERIFICA DIRETTA CL_TUNED\n');
fprintf('============================================================\n');

CL_num = ss(CL_tuned);

%% Aggiunta test
fprintf('\n============================================================\n');
fprintf(' DIAGNOSTICA POLI CL_TUNED\n');
fprintf('============================================================\n');

p_CL = pole(CL_num);

disp('Poli CL_tuned:');
disp(p_CL);

fprintf('\nParte reale dei poli:\n');
disp(real(p_CL));

fprintf('\nMassima parte reale = %.12e\n',max(real(p_CL)));

fprintf('\nNumero poli con Re >= 0:\n');
disp(sum(real(p_CL) >= 0));

fprintf('\nNumero poli con |Re| < 1e-6:\n');
disp(sum(abs(real(p_CL)) < 1e-6));

fprintf('\n============================================================\n');
fprintf(' DIAGNOSTICA INTEGRATORI\n');
fprintf('============================================================\n');

fprintf('Poli di CL_tuned vicini a zero:\n');

for k = 1:length(p_CL)
    if abs(p_CL(k)) < 1e-4
        fprintf('  pole = %.12e %+.12ei\n', ...
            real(p_CL(k)),imag(p_CL(k)));
    end
end

%% Fine aggiunta test

fprintf('Dimensione CL_tuned: %d x %d\n', ...
    size(CL_num,1),size(CL_num,2));

fprintf('Stabilità CL_tuned: %d\n',isstable(CL_num));

[g_CL_direct,w_CL_direct] = hinfnorm(CL_num);

fprintf('H-inf norm CL_tuned = %.10f\n',g_CL_direct);
fprintf('Frequenza picco      = %.10f rad/s\n',w_CL_direct);

fprintf('\nConfronto gamma:\n');
fprintf('gamma_struct         = %.10f\n',gamma_struct);
fprintf('hinfnorm(CL_tuned)   = %.10f\n',g_CL_direct);
fprintf('differenza           = %.3e\n', ...
    abs(g_CL_direct-gamma_struct));


%% ========================================================================
% 15. ESTRAZIONE DEI BLOCCHI OTTIMIZZATI
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 15. ESTRAZIONE PARAMETRI OTTIMIZZATI\n');
fprintf('============================================================\n');

Ddec_tuned = ...
    getBlockValue(CL_tuned,'Ddec');

Kalpha_tuned = ...
    getBlockValue(CL_tuned,'Kalpha');

Kbeta_tuned = ...
    getBlockValue(CL_tuned,'Kbeta');

Falpha_ss = ...
    getBlockValue(CL_tuned,'Falpha');

Fbeta_ss = ...
    getBlockValue(CL_tuned,'Fbeta');


%% ========================================================================
% 16. VISUALIZZAZIONE PARAMETRI OTTIMIZZATI
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 16. PARAMETRI OTTIMIZZATI\n');
fprintf('============================================================\n');

disp(' ');
disp('--- Kalpha ottimizzato ---');
disp(Kalpha_tuned);

disp(' ');
disp('--- Kbeta ottimizzato ---');
disp(Kbeta_tuned);

disp(' ');
disp('--- Falpha ottimizzato ---');
Falpha_tf = tf(Falpha_ss);
disp(Falpha_tf);

disp(' ');
disp('--- Fbeta ottimizzato ---');
Fbeta_tf = tf(Fbeta_ss);
disp(Fbeta_tf);

disp(' ');
disp('--- Ddec ottimizzato (scaled) ---');
disp(Ddec_tuned);


%% ========================================================================
% 17. RICOSTRUZIONE DEL CONTROLLORE SCALATO
%
% IMPORTANTE:
%
% NON usiamo minreal qui.
%
% Vogliamo prima mantenere esattamente la struttura risultante
% dall'ottimizzazione e verificare se questa riproduce CL_tuned.
%
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 17. RICOSTRUZIONE CONTROLLORE SCALATO\n');
fprintf('============================================================\n');

Calpha_scaled_raw = ...
    ss(Falpha_ss) * ss(Kalpha_tuned);

Cbeta_scaled_raw = ...
    ss(Fbeta_ss) * ss(Kbeta_tuned);

K_struct_scaled_raw = ...
    ss(Ddec_tuned) * ...
    blkdiag(Calpha_scaled_raw,Cbeta_scaled_raw);


%% ========================================================================
% 18. RICOSTRUZIONE CLOSED LOOP DAL CONTROLLORE ESTRATTO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 18. VERIFICA CONTROLLORE RICOSTRUITO\n');
fprintf('============================================================\n');

L_struct = ...
    G_scaled * K_struct_scaled_raw;

S_struct = ...
    feedback(I2,L_struct);

T_struct = ...
    feedback(L_struct,I2);

KS_struct = ...
    K_struct_scaled_raw * S_struct;


%% ========================================================================
% 19. STABILITÀ DEI TRE OGGETTI PRINCIPALI
% ========================================================================

fprintf('\n--- Stabilità ---\n');

fprintf('S  = %d\n',isstable(S_struct));
fprintf('T  = %d\n',isstable(T_struct));
fprintf('KS = %d\n',isstable(KS_struct));


%% ========================================================================
% 20. POLI
% ========================================================================

fprintf('\n--- Poli T ---\n');
disp(pole(T_struct));

fprintf('\n--- Poli KS ---\n');
disp(pole(KS_struct));

fprintf('\n--- Poli K ---\n');
disp(pole(K_struct_scaled_raw));


%% ========================================================================
% 21. H-INFINITY NORMS DEI TRE CANALI
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 21. VERIFICA H-INFINITY CONTROLLORE RICOSTRUITO\n');
fprintf('============================================================\n');

[gWS,wWS] = ...
    hinfnorm(WS*S_struct);

[gWU,wWU] = ...
    hinfnorm(WU*KS_struct);

[gWT,wWT] = ...
    hinfnorm(WT*T_struct);

CL_check = [ ...
    WS*S_struct;
    WU*KS_struct;
    WT*T_struct];

[gALL,wALL] = ...
    hinfnorm(CL_check);


fprintf('\nGamma synthesis = %.10f\n',gamma_struct);

fprintf('\nWS*S = %.10f\n',gWS);
fprintf('peak frequency = %.10f rad/s\n',wWS);

fprintf('\nWU*KS = %.10f\n',gWU);
fprintf('peak frequency = %.10f rad/s\n',wWU);

fprintf('\nWT*T = %.10f\n',gWT);
fprintf('peak frequency = %.10f rad/s\n',wWT);

fprintf('\nGamma verificato = %.10f\n',gALL);
fprintf('peak frequency   = %.10f rad/s\n',wALL);


%% ========================================================================
% 22. DC GAIN CLOSED LOOP
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 22. DC GAIN\n');
fprintf('============================================================\n');

disp('dcgain(T_struct):');
disp(dcgain(T_struct));


%% ========================================================================
% 23. CONVERSIONE ALLE UNITÀ FISICHE
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 23. CONVERSIONE IN UNITÀ FISICHE\n');
fprintf('============================================================\n');

% Relazione:
%
%       G_scaled = Dy_inv * G_nominal * Du
%
% quindi:
%
%       K_phys = Du * K_scaled * Dy_inv

K_struct = ...
    Du * K_struct_scaled_raw * Dy_inv;

% Decoupler fisico
Ddec_phys = ...
    Du * Ddec_tuned * Dy_inv;


%% ========================================================================
% 24. VERIFICA CLOSED LOOP FISICO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 24. CLOSED LOOP FISICO\n');
fprintf('============================================================\n');

L_phys = ...
    G_nominal * K_struct;

CL_phys = ...
    feedback(L_phys,I2);

fprintf('Stabilità CL fisico = %d\n',isstable(CL_phys));

fprintf('\nPoli CL fisico:\n');
disp(pole(CL_phys));

fprintf('\nDC gain CL fisico:\n');
disp(dcgain(CL_phys));


%% ========================================================================
% 25. VERSIONE MINREAL
%
% La creiamo SOLO DOPO aver fatto tutte le verifiche sopra.
%
% Questo evita che una cancellazione numerica mascheri il problema
% durante la fase di diagnostica.
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 25. VERSIONE MINREAL\n');
fprintf('============================================================\n');

K_struct_scaled = ...
    minreal(K_struct_scaled_raw,1e-7);

K_struct_minreal = ...
    minreal(K_struct,1e-7);


%% ========================================================================
% 26. INFORMAZIONI FINALI
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 26. RISULTATI FINALI\n');
fprintf('============================================================\n');

fprintf('Gamma H-inf                 = %.10f\n',gamma_struct);

fprintf('Ordine K scaled RAW         = %d\n', ...
    order(K_struct_scaled_raw));

fprintf('Ordine K scaled MINREAL     = %d\n', ...
    order(K_struct_scaled));

fprintf('Ordine K fisico MINREAL     = %d\n', ...
    order(K_struct_minreal));

fprintf('Stabilità S                 = %d\n', ...
    isstable(S_struct));

fprintf('Stabilità T                 = %d\n', ...
    isstable(T_struct));

fprintf('Stabilità KS                = %d\n', ...
    isstable(KS_struct));

fprintf('Stabilità CL fisico         = %d\n', ...
    isstable(CL_phys));


%% ========================================================================
% 27. SALVATAGGIO
% ========================================================================

fprintf('\n============================================================\n');
fprintf(' 27. SALVATAGGIO\n');
fprintf('============================================================\n');

save('HINF_controllers_struct.mat', ...
    'K_struct', ...
    'K_struct_minreal', ...
    'K_struct_scaled_raw', ...
    'K_struct_scaled', ...
    'Kalpha_tuned', ...
    'Kbeta_tuned', ...
    'Falpha_ss', ...
    'Fbeta_ss', ...
    'Falpha_tf', ...
    'Fbeta_tf', ...
    'Ddec_tuned', ...
    'Ddec_phys', ...
    'gamma_struct', ...
    'info_struct');

fprintf('\nFile salvato: HINF_controllers_struct.mat\n');


%% ========================================================================
% FINE
% ========================================================================

disp(' ');
disp('============================================================');
disp('      SINTESI HINFSTRUCT COMPLETATA');
disp('============================================================');