%% ========================================================================
% SINTESI H-INFINITO MU-SYNTHESIS (D-K ITERATION) - ELICOTTERO 2-DOF
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' MU-SYNTHESIS - IMPOSTAZIONE PESI E IMPIANTO\n');
fprintf('============================================================\n');

s = tf('s');

%% ========================================================================
% 0. CARICAMENTO DEL SETUP COMUNE H-INFINITY
% ========================================================================

% Contiene:
%   G_nominal
%   G_uncertain
%   G_scaled
%   G_uncertain_scaled
%   WS, WU, WT
%   Dy, Du, Dy_inv, Du_inv

load('HINF_workspace.mat');
load('ACTUATOR_LUMPED.mat');

fprintf('Setup H-infinity e actuator lumped caricati correttamente.\n');

%% ========================================================================
% 1. IMPIANTO INCERTO PER MU-SYNTHESIS
% ========================================================================

Gmu_scaled = G_uncertain_lumped_scaled;
Gmu_scaled.InputName = {'ubar1','ubar2'};
Gmu_scaled.OutputName = {'alpha_bar','beta_bar'};
fprintf('\nImpianto per MU-synthesis:\n');
fprintf(' Gmu_scaled = G_uncertain_lumped_scaled\n');
[~,~,blk_mu] = lftdata(Gmu_scaled);

fprintf('\nLFT USATA DA MUSYN:\n');

for k = 1:numel(blk_mu)
    fprintf('%-12s -> %2d occorrenze\n', ...
        blk_mu(k).Name, ...
        blk_mu(k).Occurrences);
end

%% ========================================================================
% CONTROLLO COERENZA NOMINALE LUMPED / NON-LUMPED
% ========================================================================

Gnom_lumped    = Gmu_scaled.NominalValue;
Gnom_nonlumped = G_unc_mu_scaled.NominalValue;

fprintf('\n============================================================\n');
fprintf('CONTROLLO COERENZA NOMINALE LUMPED / NON-LUMPED\n');
fprintf('============================================================\n');

nominalDifference = minreal( ...
    Gnom_lumped - Gnom_nonlumped, ...
    1e-8);

fprintf('||Gnom_lumped - Gnom_nonlumped||inf = %.6e\n', ...
    norm(nominalDifference,inf));

fprintf('||Gnom_lumped||inf                  = %.6e\n', ...
    norm(Gnom_lumped,inf));

fprintf('||Gnom_nonlumped||inf               = %.6e\n', ...
    norm(Gnom_nonlumped,inf));


Gphys_lumped    = G_uncertain_lumped.NominalValue;
Gphys_nonlumped = G_unc_mu.NominalValue;

physicalDifference = minreal( ...
    Gphys_lumped - Gphys_nonlumped, ...
    1e-8);

fprintf('\nConfronto modelli fisici:\n');

fprintf('||Gphys_lumped - Gphys_nonlumped||inf = %.6e\n', ...
    norm(physicalDifference,inf));

fprintf('||Gphys_lumped||inf                   = %.6e\n', ...
    norm(Gphys_lumped,inf));

fprintf('||Gphys_nonlumped||inf                = %.6e\n', ...
    norm(Gphys_nonlumped,inf));


%% ========================================================================
% 3. PESI FREQUENZIALI
% ========================================================================

% USIAMO DIRETTAMENTE I PESI DEL SETUP COMUNE:
%
%   WS
%   WU
%   WT
%
% già definiti in HINF_workspace.mat.
%
% In questo modo mixsyn, hinfsyn, hinfstruct e musyn
% lavorano con la stessa impostazione di progetto.

WS.InputName  = {'e1','e2'};
WS.OutputName = {'zS1','zS2'};

WU.InputName  = {'ubar1','ubar2'};
WU.OutputName = {'zU1','zU2'};

WT.InputName  = {'alpha_bar','beta_bar'};
WT.OutputName = {'zT1','zT2'};

%% ========================================================================
% 4. NODI DI ERRORE
% ========================================================================

% Anche i riferimenti devono essere intesi nelle variabili normalizzate.

SumE1 = sumblk('e1 = r1 - alpha_bar');
SumE2 = sumblk('e2 = r2 - beta_bar');

%% ========================================================================
% 5. IMPIANTO GENERALIZZATO P_mu
% ========================================================================

% Ingressi esogeni:
%   r1, r2      -> riferimenti normalizzati
%
% Segnali per il controllore:
%   e1, e2
%
% Comandi del controllore:
%   ubar1, ubar2
%
% Uscite pesate:
%   zS1,zS2,zU1,zU2,zT1,zT2

P_mu = connect( ...
    Gmu_scaled, ...
    WS, ...
    WU, ...
    WT, ...
    SumE1, ...
    SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2','e1','e2'});

%% ========================================================================
% 6. D-K ITERATION - MUSYN
% ========================================================================

nmeas = 2;
ncont = 2;

optsMU = musynOptions( ...
    'Display','short', ...
    'MixedMU','on', ...
    'FullDG',true, ...
    'TargetPerf',0.99, ...
    'MaxIter',40, ...
    'TolPerf',0);

fprintf('\n============================================================\n');
fprintf('AVVIO D-K ITERATION\n');
fprintf('Obiettivo: Robust Performance < 1\n');
fprintf('============================================================\n');

[K_mu_scaled, CLperf_mu, info_mu] = ...
    musyn(P_mu,nmeas,ncont,optsMU);

fprintf('\n============================================================\n');
fprintf('VERIFICA CONTROLLORE RAW\n');
fprintf('============================================================\n');

Kraw = ss(K_mu_scaled);

fprintf('Ordine K raw = %d\n',order(Kraw));
fprintf('K raw stabile = %d\n',isstable(Kraw));

G_nom_scaled = Gmu_scaled.NominalValue;
I2 = eye(2);

Lraw = G_nom_scaled * Kraw;
Sraw = feedback(I2,Lraw);
Traw = feedback(Lraw,I2);

fprintf('Closed-loop nominale raw stabile = %d\n',isstable(Traw));
fprintf('Max Re polo raw = %.6e\n',max(real(pole(Traw))));

%% ========================================================================
% 7. VERIFICA E RIDUZIONE DELL'ORDINE SUL MODELLO NON-LUMPED
% ========================================================================

% La sintesi MU viene effettuata sul modello lumped:
%
%     G_uncertain_lumped_scaled
%
% La validazione finale della Robust Performance viene invece effettuata
% sul modello parametrico NON-LUMPED:
%
%     G_unc_mu_scaled
%
% contenente:
%     J_alpha, l, omega_n, tau_d
%
% Pertanto anche la riduzione dell'ordine viene accettata
% solo se mantiene RP < 1 sul modello NON-LUMPED.

fprintf('\n============================================================\n');
fprintf('VERIFICA RP E RIDUZIONE SUL MODELLO NON-LUMPED\n');
fprintf('============================================================\n');

% Controllore full-order prodotto dalla mu-synthesis
K_mu_scaled_full = ss(K_mu_scaled);

N = order(K_mu_scaled_full);

fprintf('Ordine controllore full-order = %d\n',N);
fprintf('RP della sintesi sul lumped    = %.6f\n',CLperf_mu);


%% ------------------------------------------------------------------------
% 7.1 COSTRUZIONE CLOSED-LOOP NON-LUMPED FULL-ORDER
% -------------------------------------------------------------------------

Gmu_nonlumped = G_unc_mu_scaled;

Gmu_nonlumped.InputName  = {'u1','u2'};
Gmu_nonlumped.OutputName = {'y1','y2'};

Kfull = K_mu_scaled_full;
Kfull.InputName  = {'e1','e2'};
Kfull.OutputName = {'u1','u2'};

W_S = ss(WS);
W_S.InputName  = {'e1','e2'};
W_S.OutputName = {'zS1','zS2'};

W_U = ss(WU);
W_U.InputName  = {'u1','u2'};
W_U.OutputName = {'zU1','zU2'};

W_T = ss(WT);
W_T.InputName  = {'y1','y2'};
W_T.OutputName = {'zT1','zT2'};

sum1 = sumblk('e1 = r1 - y1');
sum2 = sumblk('e2 = r2 - y2');

CL_full_nonlumped = connect( ...
    Gmu_nonlumped, ...
    Kfull, ...
    W_S, ...
    W_U, ...
    W_T, ...
    sum1, ...
    sum2, ...
    {'r1','r2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2'});


%% ------------------------------------------------------------------------
% 7.2 ROBUST PERFORMANCE FULL-ORDER NON-LUMPED
% -------------------------------------------------------------------------

fprintf('\nVerifica RP full-order sul modello NON-LUMPED...\n');

optsWC = wcOptions( ...
    'Display','off', ...
    'MussvOptions','a');

[wc_full,wcu_full] = wcgain( ...
    CL_full_nonlumped, ...
    optsWC);
fprintf('\n============================================================\n');
fprintf('WORST-CASE UNCERTAINTIES\n');
fprintf('============================================================\n');

fprintf('Worst-case J_alpha = %.12g\n', wcu_full.J_alpha);
fprintf('Worst-case l       = %.12g\n', wcu_full.l);
fprintf('Worst-case omega_n = %.12g\n', wcu_full.omega_n);
fprintf('Worst-case tau_d   = %.12g\n', wcu_full.tau_d);

RP_full_nonlumped = wc_full.UpperBound;

fprintf('RP full-order NON-LUMPED lower = %.6f\n', ...
    wc_full.LowerBound);

fprintf('RP full-order NON-LUMPED upper = %.6f\n', ...
    wc_full.UpperBound);

fprintf('Frequenza critica = %.6f rad/s\n', ...
    wc_full.CriticalFrequency);



fprintf('RP full-order NON-LUMPED = %.6f\n', ...
    RP_full_nonlumped);

if RP_full_nonlumped < 1

    fprintf('FULL-ORDER: RP < 1 [PASS]\n');

else

    fprintf('FULL-ORDER: RP >= 1 [FAIL]\n');

    error(['Il controllore full-order ottenuto dalla mu-synthesis ' ...
           'NON soddisfa RP < 1 sul modello non-lumped. ' ...
           'La riduzione d''ordine non deve essere eseguita. ' ...
           'Occorre modificare la sintesi.']);

end

fprintf('\n============================================================\n');
fprintf('RIDUZIONE OTTIMIZZATA DELL''ORDINE\n');
fprintf('============================================================\n');

K_mu_scaled_reduced = [];
RP_reduced = NaN;
order_reduced = NaN;

for k = 8:N

    Ktest = ss(Kred(:,:,k));

    Ktest.InputName  = {'e1','e2'};
    Ktest.OutputName = {'u1','u2'};

    % ---------------------------------------------------------
    % Closed-loop incerto NON-LUMPED
    % ---------------------------------------------------------
    CLtest = connect( ...
        Gmu_nonlumped, ...
        Ktest, ...
        W_S, ...
        W_U, ...
        W_T, ...
        sum1, ...
        sum2, ...
        {'r1','r2'}, ...
        {'zS1','zS2','zU1','zU2','zT1','zT2'});

    % ---------------------------------------------------------
    % 1) Controllo NOMINALE
    % ---------------------------------------------------------
    CLnom_test = nominal(CLtest);

    if ~isstable(CLnom_test)

        fprintf('Ordine %3d --> NP = Inf --> skip\n',k);
        continue;

    end

    NP_test = hinfnorm( ...
        minreal(CLnom_test,1e-5));

    % Se NP >= 1, RP non può essere < 1
    if NP_test >= 1

        fprintf('Ordine %3d --> NP = %.6f --> skip WCgain\n', ...
            k,NP_test);

        continue;

    end

    % ---------------------------------------------------------
    % 2) Solo se NP < 1 facciamo wcgain
    % ---------------------------------------------------------
    [wc_test,~] = wcgain( ...
        CLtest, ...
        optsWC);

    RP_test = wc_test.UpperBound;

    fprintf('Ordine %3d --> NP = %.6f | RP = %.6f\n', ...
        k,NP_test,RP_test);

    % Primo ordine che soddisfa RP < 1
    if RP_test < 1

        K_mu_scaled_reduced = Ktest;
        RP_reduced = RP_test;
        order_reduced = k;

        fprintf('\nPASS trovato!\n');
        fprintf('Ordine minimo trovato = %d\n',order_reduced);
        fprintf('RP = %.6f\n',RP_reduced);

        break;
    end
end

%% ------------------------------------------------------------------------
% 7.4 SELEZIONE CONTROLLore FINALE
% -------------------------------------------------------------------------

if isempty(K_mu_scaled_reduced)

    fprintf('\nNESSUNA RIDUZIONE SODDISFA RP < 1\n');
    fprintf('Viene mantenuto il controllore full-order.\n');

    K_mu_scaled = K_mu_scaled_full;

    RP_final = RP_full_nonlumped;
    order_final = N;

else

    fprintf('\nRIDUZIONE RIUSCITA SUL MODELLO NON-LUMPED.\n');

    K_mu_scaled = K_mu_scaled_reduced;

    RP_final = RP_reduced;
    order_final = order_reduced;

    fprintf('Ordine controllore finale = %d\n', ...
        order_final);

    fprintf('RP finale NON-LUMPED       = %.6f\n', ...
        RP_final);
end
%% ========================================================================
% 8. DESCALATURA DEL CONTROLLORE
% ========================================================================

% La relazione di normalizzazione è:
%
%       u = Du * u_bar
%       e_bar = Dy^-1 * e
%
% quindi:
%
%       K_phys = Du * K_scaled * Dy^-1

K_mu = Du * K_mu_scaled * Dy_inv;
% Nomi per eventuale utilizzo in connect / Simulink

K_mu.InputName  = {'e1','e2'};
K_mu.OutputName = {'u1','u2'};

fprintf('Ordine controllore finale    = %d\n', ...
    order(K_mu));

fprintf('\nControllore riportato alle unita'' fisiche.\n');

%% ========================================================================
% 9. VERIFICA DEL RISULTATO MU
% ========================================================================

fprintf('\n============================================================\n');
fprintf('VERIFICA RISULTATO MU\n');
fprintf('============================================================\n');

fprintf('Ordine controllore finale = %d\n', ...
    order(K_mu_scaled));

fprintf('Robust Performance finale = %.6f\n', ...
    RP_final);

if RP_final < 1

    fprintf('\nSUCCESSO!\n');
    fprintf(['Il controllore finale soddisfa la specifica ', ...
             'Robust Performance < 1.\n']);

else

    fprintf('\nATTENZIONE!\n');
    fprintf(['Robust Performance >= 1. ', ...
             'Le specifiche robuste non sono soddisfatte.\n']);

end

%% ========================================================================
% 10. VERIFICA NOMINALE DELLA PRESTAZIONE
% ========================================================================

% IMPORTANTE:
% la verifica nominale delle prestazioni pesate viene fatta
% sul problema SCALATO, quindi usiamo:
%
%   Gmu_scaled
%   K_mu_scaled
%
% e non Gmu + K_mu direttamente.

G_nom_nonlumped = G_unc_mu_scaled.NominalValue;

Lnom_nonlumped = G_nom_nonlumped * K_mu_scaled;
Snom_nonlumped = feedback(I2,Lnom_nonlumped);
Tnom_nonlumped = feedback(Lnom_nonlumped,I2);
KSnom_nonlumped = K_mu_scaled * Snom_nonlumped;

CLnom_nonlumped = [ ...
    WS*Snom_nonlumped;
    WU*KSnom_nonlumped;
    WT*Tnom_nonlumped];

gammaNom_nonlumped = hinfnorm( ...
    minreal(CLnom_nonlumped,1e-5));

fprintf('Gamma nominale NON-LUMPED = %.6f\n', ...
    gammaNom_nonlumped);

%% ========================================================================
% 10.1 NOMINAL PERFORMANCE ESATTA DEL PROBLEMA MU
% ========================================================================

CLnom_MU = lft( ...
    P_mu.NominalValue, ...
    K_mu_scaled);

CLnom_MU = minreal(CLnom_MU,1e-5);

gammaNom_MU = hinfnorm(CLnom_MU);

fprintf('\n============================================================\n');
fprintf('NOMINAL PERFORMANCE - PROBLEMA MU\n');
fprintf('============================================================\n');

fprintf('NP MU = %.6f\n',gammaNom_MU);

if gammaNom_MU < 1
    fprintf('NP MU < 1 [PASS]\n');
else
    fprintf('NP MU >= 1 [FAIL]\n');
end

%% ========================================================================
% 11. VERIFICA NOMINAL STABILITY (NS)
% ========================================================================

fprintf('\n============================================================\n');
fprintf('VERIFICA NOMINAL STABILITY (NS)\n');
fprintf('============================================================\n');

% Verifica nominale sul modello NON-LUMPED
poles_nominal = pole(Tnom_nonlumped);
maxRePole = max(real(poles_nominal));

fprintf('Max Re(polo closed-loop nominale NON-LUMPED) = %.6e\n', ...
    maxRePole);

if maxRePole < 0

    fprintf('NS = OK: sistema nominalmente stabile.\n');

else

    fprintf('NS = ATTENZIONE: sistema nominalmente instabile.\n');

end

%% ========================================================================
% 12. VERIFICA ROBUST STABILITY (RS)
% ========================================================================

% Verifica sul modello parametrico non-lumped.
%
% Il controllore e' stato sintetizzato sul modello lumped:
%
%     G_uncertain_lumped
%
% ma la verifica viene effettuata sul modello parametrico originale
% contenente le quattro incertezze:
%
%     J_alpha, l, omega_n, tau_d
%
Gmu_phys = G_unc_mu;

Gmu_phys.InputName  = {'u1','u2'};
Gmu_phys.OutputName = {'alpha','beta'};

SumE1_phys = sumblk('e1 = r1 - alpha');
SumE2_phys = sumblk('e2 = r2 - beta');

CL_unc_mu = connect( ...
    Gmu_phys, ...
    K_mu, ...
    SumE1_phys, ...
    SumE2_phys, ...
    {'r1','r2'}, ...
    {'alpha','beta'});

[stabmarg,~,info_RS] = robuststab(CL_unc_mu);

fprintf('Robust Stability Lower Bound = %.6f\n', ...
    stabmarg.LowerBound);

fprintf('Robust Stability Upper Bound = %.6f\n', ...
    stabmarg.UpperBound);

if stabmarg.LowerBound > 1

    fprintf('RS = OK: il sistema e'' robustamente stabile.\n');

else

    fprintf('RS = ATTENZIONE: robust stability non garantita.\n');

end

%% ========================================================================
% 13. SALVATAGGIO
% ========================================================================

save('MU_controller.mat', ...
    'K_mu', ...
    'K_mu_scaled', ...
    'K_mu_scaled_full', ...
    'CLperf_mu', ...
    'RP_full_nonlumped', ...
    'RP_final', ...
    'gammaNom_nonlumped', ...
    'order_final', ...
    'stabmarg', ...
    'info_mu', ...
    'info_RS', ...
    'P_mu');

fprintf('\n============================================================\n');
fprintf('SALVATAGGIO COMPLETATO\n');
fprintf('============================================================\n');

fprintf('File "MU_controller.mat" creato correttamente.\n');