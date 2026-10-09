%% ========================================================================
% Hinf_sintesi.m
%
% Sintesi H-infinity per l'elicottero 2DoF
%
% Confronto tra:
%
%   1) mixsyn
%   2) hinfsyn
%   3) hinfstruct con struttura:
%
%          PID_alpha -> Lead/Lag_alpha --+
%                                        |
%          PID_beta  -> Lead/Lag_beta  --+-> Ddec -> u
%
% Il controllore strutturato e':
%
%   K(s) = Ddec *
%          blkdiag(Falpha*Kalpha, Fbeta*Kbeta)
%
% La sintesi viene effettuata sul plant normalizzato G_scaled.
% Il controllore viene successivamente riportato nelle coordinate fisiche.
%
% ========================================================================

close all;
clc;

%% ========================================================================
% 0.2 PARAMETRI GENERALI
% ========================================================================

I2 = eye(2);

nmeas = 2;
ncont = 2;

% Variabile di Laplace.
% Necessaria per costruire i compensatori iniziali.
s = tf('s');

% Tolleranza per minreal.
tol = 1e-7;

% Opzioni per mixsyn / hinfsyn.
optsHinf = hinfsynOptions( ...
    'Display', 'on');


%% ========================================================================
% 1. MIXSYN
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('1. MIXSYN\n');
fprintf('============================================================\n');

[K_mix_scaled, ...
 CL_mix, ...
 gamma_mix, ...
 info_mix] = ...
    mixsyn( ...
        G_scaled, ...
        WS, ...
        WU, ...
        WT, ...
        optsHinf);

if isempty(K_mix_scaled) || ~isfinite(gamma_mix)

    error('La sintesi mixsyn non ha prodotto un controllore valido.');

end

K_mix_scaled = ...
    minreal( ...
        ss(K_mix_scaled), ...
        tol);


% Conversione nelle coordinate fisiche:
%
%   y = Dy*y_scaled
%   u = Du*u_scaled
%
% quindi:
%
%   K_phys = Du*K_scaled*Dy_inv

K_mix = ...
    minreal( ...
        Du * ...
        K_mix_scaled * ...
        Dy_inv, ...
        tol);


fprintf('Gamma mixsyn       = %.8f\n', gamma_mix);
fprintf('Ordine K_mixsyn    = %d\n', order(K_mix_scaled));


%% ========================================================================
% 2. HINFSYN
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('2. HINFSYN\n');
fprintf('============================================================\n');

[K_hinfsyn_scaled, ...
 CL_hinfsyn, ...
 gamma_hinfsyn, ...
 info_hinfsyn] = ...
    hinfsyn( ...
        P_mix, ...
        nmeas, ...
        ncont, ...
        optsHinf);

if isempty(K_hinfsyn_scaled) || ~isfinite(gamma_hinfsyn)

    error('La sintesi hinfsyn non ha prodotto un controllore valido.');

end

K_hinfsyn_scaled = ...
    minreal( ...
        ss(K_hinfsyn_scaled), ...
        tol);


K_hinfsyn = ...
    minreal( ...
        Du * ...
        K_hinfsyn_scaled * ...
        Dy_inv, ...
        tol);


fprintf('Gamma hinfsyn      = %.8f\n', gamma_hinfsyn);
fprintf('Ordine K_hinfsyn   = %d\n', order(K_hinfsyn_scaled));

fprintf( ...
    'Differenza relativa gamma = %.6e\n', ...
    abs(gamma_mix - gamma_hinfsyn) / max(gamma_mix, eps));


%% ========================================================================
% 3. HINFSTRUCT
%
% PID + LEAD/LAG + DISACCOPPIATORE
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('3. HINFSTRUCT - PID + LEAD/LAG + DDEC\n');
fprintf('============================================================\n');


%% ========================================================================
% 3.1 NOMI DEI SEGNALI
% ========================================================================

G_scaled.InputName = {
    'u1'
    'u2'
};

G_scaled.OutputName = {
    'y1'
    'y2'
};

WS.InputName = {
    'e1'
    'e2'
};

WS.OutputName = {
    'zS1'
    'zS2'
};

WU.InputName = {
    'u1'
    'u2'
};

WU.OutputName = {
    'zU1'
    'zU2'
};

WT.InputName = {
    'y1'
    'y2'
};

WT.OutputName = {
    'zT1'
    'zT2'
};


%% ========================================================================
% 3.2 PID TUNABLE
% ========================================================================
%
% Il PID e' inizializzato usando i valori fisici gia' utilizzati
% nel progetto e convertendoli nelle coordinate normalizzate.
%
% La conversione vale:
%
%   K_scaled = (scale_angle / scale_force) * K_phys
%
% per ciascun canale.

alphaPhysToScaled = ...
    scale_alpha / scale_F1;

betaPhysToScaled = ...
    scale_beta / scale_F2;


%% PID alpha

Kp_alpha_phys_0 = -0.2147;
Ki_alpha_phys_0 =  1.4066;

Kalpha_tunable = ...
    tunablePID( ...
        'Kalpha', ...
        'PID');

Kalpha_tunable.Kp.Value = ...
    alphaPhysToScaled * ...
    Kp_alpha_phys_0;

Kalpha_tunable.Ki.Value = ...
    alphaPhysToScaled * ...
    Ki_alpha_phys_0;

Kalpha_tunable.Kd.Value = ...
    0;

Kalpha_tunable.Tf.Value = ...
    0.02;


%% PID beta

Kp_beta_phys_0 = 0.2149;
Ki_beta_phys_0 = 0.0444;

Kbeta_tunable = ...
    tunablePID( ...
        'Kbeta', ...
        'PID');

Kbeta_tunable.Kp.Value = ...
    betaPhysToScaled * ...
    Kp_beta_phys_0;

Kbeta_tunable.Ki.Value = ...
    betaPhysToScaled * ...
    Ki_beta_phys_0;

Kbeta_tunable.Kd.Value = ...
    0;

Kbeta_tunable.Tf.Value = ...
    0.02;


%% Vincoli PID

Kalpha_tunable.Ki.Minimum = 1e-10;
Kbeta_tunable.Ki.Minimum  = 1e-10;

Kalpha_tunable.Tf.Minimum = 1e-4;
Kbeta_tunable.Tf.Minimum  = 1e-4;


%% ========================================================================
% 3.3 COMPENSATORI LEAD/LAG TUNABLE
% ========================================================================
%
% I pesi di tracking hanno:
%
%   wb_alpha = 3.8 rad/s
%   wb_beta  = 3.0 rad/s
%
% mentre WT interviene principalmente a frequenze piu' alte:
%
%   wt_alpha = 22 rad/s
%   wt_beta  = 18 rad/s
%
% Per questo non usiamo direttamente i valori arbitrari dell'esempio
% dei colleghi.
%
% L'inizializzazione e' scelta in prossimita' della banda di tracking,
% lasciando comunque liberi i coefficienti durante hinfstruct.
%
% Forma:
%
%              1 + s/wz
% F(s) = ---------------------
%              1 + s/wp
%
% Inizialmente:
%
% alpha: wz = 3.8, wp = 12
% beta : wz = 3.0, wp = 10
%
%

wz_alpha_0 = 3.8;
wp_alpha_0 = 12.0;

wz_beta_0 = 3.0;
wp_beta_0 = 10.0;


Falpha0 = ...
    (1 + s/wz_alpha_0) / ...
    (1 + s/wp_alpha_0);

Fbeta0 = ...
    (1 + s/wz_beta_0) / ...
    (1 + s/wp_beta_0);


Falpha_tunable = ...
    tunableTF( ...
        'Falpha', ...
        tf(Falpha0));

Fbeta_tunable = ...
    tunableTF( ...
        'Fbeta', ...
        tf(Fbeta0));


%% Vincoli di stabilita' dei compensatori
%
% tunableTF ha denominatore normalizzato con primo coefficiente = 1.
% Il termine costante deve rimanere positivo per mantenere il polo
% reale nel semipiano sinistro.

Falpha_tunable.Denominator.Minimum(end) = 1e-4;
Fbeta_tunable.Denominator.Minimum(end)  = 1e-4;


%% ========================================================================
% 3.4 PID + LEAD/LAG
% ========================================================================

Calpha_tunable = ...
    Falpha_tunable * ...
    Kalpha_tunable;

Cbeta_tunable = ...
    Fbeta_tunable * ...
    Kbeta_tunable;


%% ========================================================================
% 3.5 STRUTTURA DIAGONALE
% ========================================================================
%
%       [ Falpha*Kalpha       0       ]
% Kd =  [                           ]
%       [      0        Fbeta*Kbeta  ]

Kdiag_tunable = ...
    blkdiag( ...
        Calpha_tunable, ...
        Cbeta_tunable);


%% ========================================================================
% 3.6 DISACCOPPIATORE STATICO
% ========================================================================

Gdc_scaled = dcgain(G_scaled);

if ~all(isfinite(Gdc_scaled(:))) || ...
        rcond(Gdc_scaled) < 1e-8

    warning( ...
        ['dcgain(G_scaled) non valido o mal condizionato. ' ...
         'Uso la risposta a bassa frequenza a 0.1 rad/s.']);

    w_ref = 0.1;

    Gdc_scaled = ...
        real(freqresp(G_scaled, w_ref));

end


Ddec0 = ...
    Gdc_scaled \ eye(2);


Ddec_tunable = ...
    tunableGain( ...
        'Ddec', ...
        2, ...
        2);

Ddec_tunable.Gain.Value = ...
    Ddec0;


%% ========================================================================
% 3.7 CONTROLLORE STRUTTURATO COMPLETO
% ========================================================================
%
%       e
%       |
%       +--> PID_alpha --> Falpha --+
%                                   |
%       +--> PID_beta  --> Fbeta ---+--> Ddec --> u

K_pidcomp_tunable = ...
    Ddec_tunable * ...
    Kdiag_tunable;

K_pidcomp_tunable.InputName = {
    'e1'
    'e2'
};

K_pidcomp_tunable.OutputName = {
    'u1'
    'u2'
};


%% ========================================================================
% 3.8 CLOSED LOOP TUNABLE
% ========================================================================

sum1 = ...
    sumblk('e1 = r1 - y1');

sum2 = ...
    sumblk('e2 = r2 - y2');


CL_pidcomp_tunable = ...
    connect( ...
        G_scaled, ...
        WS, ...
        WU, ...
        WT, ...
        K_pidcomp_tunable, ...
        sum1, ...
        sum2, ...
        {'r1','r2'}, ...
        {'zS1','zS2', ...
         'zU1','zU2', ...
         'zT1','zT2'});


%% ========================================================================
% 3.9 HINFSTRUCT
% ========================================================================

rng(10);

optsHS = hinfstructOptions( ...
    'Display',      'final', ...
    'RandomStart',  10, ...
    'MaxIter',      500, ...
    'TargetGain',   0, ...
    'MinDecay',     0.02);

[CL_pidcomp_tuned, ...
 gamma_pidcomp, ...
 info_pidcomp] = ...
    hinfstruct( ...
        CL_pidcomp_tunable, ...
        optsHS);


if isempty(CL_pidcomp_tuned) || ...
        ~isfinite(gamma_pidcomp)

    error('hinfstruct non ha prodotto un risultato valido.');

end

%% Verifica delle partenze

fprintf('\n');
fprintf('Numero di ottimizzazioni eseguite: %d\n', ...
    numel(info_pidcomp));

fprintf('Miglior gamma di sintesi: %.8f\n', ...
    gamma_pidcomp);

fprintf('\nGamma delle singole ottimizzazioni, ordinati:\n');
disp(sort([info_pidcomp.Objective]).');


%% ========================================================================
% 3.10 ESTRAZIONE DEI PARAMETRI OTTIMIZZATI
% ========================================================================

Kalpha_tuned_scaled = ...
    getBlockValue( ...
        CL_pidcomp_tuned, ...
        'Kalpha');

Kbeta_tuned_scaled = ...
    getBlockValue( ...
        CL_pidcomp_tuned, ...
        'Kbeta');


Falpha_tuned = ...
    getBlockValue( ...
        CL_pidcomp_tuned, ...
        'Falpha');

Fbeta_tuned = ...
    getBlockValue( ...
        CL_pidcomp_tuned, ...
        'Fbeta');


Ddec_tuned = ...
    getBlockValue( ...
        CL_pidcomp_tuned, ...
        'Ddec');


Ddec_scaled = ...
    dcgain(Ddec_tuned);


%% ========================================================================
% 3.11 CONTROLLORE NUMERICO SCALATO
% ========================================================================

Calpha_scaled = ...
    minreal( ...
        ss(Falpha_tuned) * ...
        ss(Kalpha_tuned_scaled), ...
        tol);

Cbeta_scaled = ...
    minreal( ...
        ss(Fbeta_tuned) * ...
        ss(Kbeta_tuned_scaled), ...
        tol);


Kdiag_scaled = ...
    blkdiag( ...
        Calpha_scaled, ...
        Cbeta_scaled);


K_pidcomp_scaled = ...
    minreal( ...
        ss(Ddec_tuned) * ...
        Kdiag_scaled, ...
        tol);


%% ========================================================================
% 3.12 CONVERSIONE IN COORDINATE FISICHE
% ========================================================================

K_pidcomp = ...
    minreal( ...
        Du * ...
        K_pidcomp_scaled * ...
        Dy_inv, ...
        tol);


%% ========================================================================
% 3.13 PID FISICI
% ========================================================================

alphaScaledToPhys = ...
    scale_F1 / scale_alpha;

betaScaledToPhys = ...
    scale_F2 / scale_beta;


Kalpha_tuned = ...
    pid( ...
        alphaScaledToPhys * Kalpha_tuned_scaled.Kp, ...
        alphaScaledToPhys * Kalpha_tuned_scaled.Ki, ...
        alphaScaledToPhys * Kalpha_tuned_scaled.Kd, ...
        Kalpha_tuned_scaled.Tf);


Kbeta_tuned = ...
    pid( ...
        betaScaledToPhys * Kbeta_tuned_scaled.Kp, ...
        betaScaledToPhys * Kbeta_tuned_scaled.Ki, ...
        betaScaledToPhys * Kbeta_tuned_scaled.Kd, ...
        Kbeta_tuned_scaled.Tf);


%% ========================================================================
% 3.14 PARAMETRI PID PER SIMULINK
% ========================================================================

Kp_alpha = Kalpha_tuned.Kp;
Ki_alpha = Kalpha_tuned.Ki;
Kd_alpha = Kalpha_tuned.Kd;
Tf_alpha = Kalpha_tuned.Tf;

Kp_beta = Kbeta_tuned.Kp;
Ki_beta = Kbeta_tuned.Ki;
Kd_beta = Kbeta_tuned.Kd;
Tf_beta = Kbeta_tuned.Tf;


%% ========================================================================
% 3.15 COMPENSATORI PER SIMULINK
% ========================================================================
%
% Forma:
%
%       b1*s + b0
% F = --------------
%          s + a0
%
% Il primo coefficiente del denominatore e' quindi 1.

[numFalpha, denFalpha] = ...
    tfdata( ...
        tf(Falpha_tuned), ...
        'v');

[numFbeta, denFbeta] = ...
    tfdata( ...
        tf(Fbeta_tuned), ...
        'v');


%% ========================================================================
% 3.16 DISACCOPPIATORE FISICO
% ========================================================================

Ddec_phys = ...
    Du * ...
    Ddec_scaled * ...
    inv(Du);


%% ========================================================================
% 3.17 VERIFICA DELLA SINTESI
% ========================================================================

L_pidcomp = ...
    G_scaled * ...
    K_pidcomp_scaled;

S_pidcomp = ...
    feedback( ...
        I2, ...
        L_pidcomp);

T_pidcomp = ...
    feedback( ...
        L_pidcomp, ...
        I2);

KS_pidcomp = ...
    K_pidcomp_scaled * ...
    S_pidcomp;


%% Norme singole

[gWS_pidcomp, wWS_pidcomp] = ...
    hinfnorm( ...
        minreal( ...
            WS * S_pidcomp, ...
            tol));


[gWU_pidcomp, wWU_pidcomp] = ...
    hinfnorm( ...
        minreal( ...
            WU * KS_pidcomp, ...
            tol));


[gWT_pidcomp, wWT_pidcomp] = ...
    hinfnorm( ...
        minreal( ...
            WT * T_pidcomp, ...
            tol));


%% Norma complessiva

CLcheck = [
    WS * S_pidcomp
    WU * KS_pidcomp
    WT * T_pidcomp
];


[gALL_pidcomp, wALL_pidcomp] = ...
    hinfnorm( ...
        minreal( ...
            CLcheck, ...
            tol));


%% ========================================================================
% 3.18 VERIFICA CONVERSIONE FISICA
% ========================================================================
%
% Questa verifica e' valida perche' il controllore completo e'
%
%   Kphys = Du*Kscaled*Dy_inv.
%
% La realizzazione fisica separata e':
%
%   Ddec_phys *
%   blkdiag(PID_alpha_phys, PID_beta_phys)
%
% con i compensatori dinamici invariati.




%% ========================================================================
% 3.19 REPORT
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('RISULTATO HINFSTRUCT\n');
fprintf('============================================================\n');

fprintf('Gamma synthesis = %.8f\n', gamma_pidcomp);
fprintf('Gamma verified  = %.8f\n', gALL_pidcomp);

fprintf('\n');
fprintf('Contributi:\n');

fprintf('  ||WS*S||inf  = %.8f\n', gWS_pidcomp);
fprintf('  ||WU*KS||inf = %.8f\n', gWU_pidcomp);
fprintf('  ||WT*T||inf  = %.8f\n', gWT_pidcomp);

fprintf('\n');
fprintf('Massima parte reale dei poli = %.6e\n', ...
    max(real(pole(T_pidcomp))));

%% PID alpha

fprintf('\n------------------------------------------------------------\n');
fprintf('PID ALPHA - FISICO\n');
fprintf('------------------------------------------------------------\n');

disp(Kalpha_tuned);


%% PID beta

fprintf('\n------------------------------------------------------------\n');
fprintf('PID BETA - FISICO\n');
fprintf('------------------------------------------------------------\n');

disp(Kbeta_tuned);


%% Compensatore alpha

fprintf('\n------------------------------------------------------------\n');
fprintf('LEAD/LAG ALPHA\n');
fprintf('------------------------------------------------------------\n');

disp(tf(Falpha_tuned));

fprintf('Numeratore:\n');
disp(numFalpha);

fprintf('Denominatore:\n');
disp(denFalpha);


%% Compensatore beta

fprintf('\n------------------------------------------------------------\n');
fprintf('LEAD/LAG BETA\n');
fprintf('------------------------------------------------------------\n');

disp(tf(Fbeta_tuned));

fprintf('Numeratore:\n');
disp(numFbeta);

fprintf('Denominatore:\n');
disp(denFbeta);


%% Disaccoppiatore

fprintf('\n------------------------------------------------------------\n');
fprintf('DDEC SCALATO\n');
fprintf('------------------------------------------------------------\n');

disp(Ddec_scaled);


fprintf('\n------------------------------------------------------------\n');
fprintf('DDEC FISICO\n');
fprintf('------------------------------------------------------------\n');

disp(Ddec_phys);


%% ========================================================================
% 4. CONTROLLORI FINALI
% ========================================================================

controllersScaled = {
    K_mix_scaled
    K_hinfsyn_scaled
    K_pidcomp_scaled
};


controllersPhysical = {
    K_mix
    K_hinfsyn
    K_pidcomp
};


controllerNames = {
    'mixsyn'
    'hinfsyn'
    'PID + Lead/Lag + Ddec'
};


nControllers = ...
    numel(controllersScaled);


%% ========================================================================
% 5. SALVATAGGIO
% ========================================================================

save( ...
    'HINF_controllers.mat', ...
    'K_mix_scaled', ...
    'K_hinfsyn_scaled', ...
    'K_pidcomp_scaled', ...
    'K_mix', ...
    'K_hinfsyn', ...
    'K_pidcomp', ...
    'Kalpha_tuned_scaled', ...
    'Kbeta_tuned_scaled', ...
    'Kalpha_tuned', ...
    'Kbeta_tuned', ...
    'Falpha_tuned', ...
    'Fbeta_tuned', ...
    'numFalpha', ...
    'denFalpha', ...
    'numFbeta', ...
    'denFbeta', ...
    'gamma_mix', ...
    'gamma_hinfsyn', ...
    'gamma_pidcomp', ...
    'info_mix', ...
    'info_hinfsyn', ...
    'info_pidcomp', ...
    'nControllers', ...
    'controllersScaled', ...
    'controllersPhysical', ...
    'controllerNames', ...
    'Ddec_scaled', ...
    'Ddec_phys', ...
    'Kp_alpha', ...
    'Ki_alpha', ...
    'Kd_alpha', ...
    'Tf_alpha', ...
    'Kp_beta', ...
    'Ki_beta', ...
    'Kd_beta', ...
    'Tf_beta', ...
    'gWS_pidcomp', ...
    'gWU_pidcomp', ...
    'gWT_pidcomp', ...
    'gALL_pidcomp');


fprintf('\n');
fprintf('============================================================\n');
fprintf('HINF_controllers.mat creato correttamente.\n');
fprintf('============================================================\n');
