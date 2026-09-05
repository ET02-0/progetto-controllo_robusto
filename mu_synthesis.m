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

% Uniamo la cinematica incerta con la dinamica degli attuatori incerta

Gmu_scaled = G_uncertain_lumped_scaled;
Gmu_scaled.InputName  = {'ubar1','ubar2'};
Gmu_scaled.OutputName = {'alpha_bar','beta_bar'};



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
    'Display','full', ...
    'MixedMU','off', ...
    'MaxIter',2);

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
% 7. RIDUZIONE DELL'ORDINE DEL CONTROLLORE SCALATO
% ========================================================================

K_mu_scaled = ss(K_mu_scaled);

fprintf('\n============================================================\n');
fprintf('RISULTATI SINTESI MU - CONTROLLORE SCALATO\n');
fprintf('============================================================\n');

fprintf('Performance robusta (mu) = %.6f\n',CLperf_mu);

fprintf('Ordine controllore originario = %d\n', ...
    order(K_mu_scaled));

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

fprintf('Ordine controllore ridotto    = %d\n', ...
    order(K_mu));

fprintf('\nControllore riportato alle unita'' fisiche.\n');

%% ========================================================================
% 9. VERIFICA DEL RISULTATO MU
% ========================================================================

if CLperf_mu < 1

    fprintf('\nSUCCESSO!\n');
    fprintf(['La sintesi mu ha ottenuto Robust Performance < 1 ', ...
             'sul problema generalizzato scalato.\n']);

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

G_nom_scaled = Gmu_scaled.NominalValue;

I2 = eye(2);

Lmu_scaled = ...
    G_nom_scaled * K_mu_scaled;

Smu_scaled = ...
    feedback(I2,Lmu_scaled);

Tmu_scaled = ...
    feedback(Lmu_scaled,I2);

KSmu_scaled = ...
    K_mu_scaled * Smu_scaled;

CLnom_mu = [ ...
    WS * Smu_scaled;
    WU * KSmu_scaled;
    WT * Tmu_scaled];

gammaNom_mu = hinfnorm( ...
    minreal(CLnom_mu,1e-5));

fprintf('\n============================================================\n');
fprintf('VERIFICA NOMINALE (NP)\n');
fprintf('============================================================\n');

fprintf('Gamma nominale = %.6f\n',gammaNom_mu);

%% ========================================================================
% 11. VERIFICA NOMINAL STABILITY (NS)
% ========================================================================

fprintf('\n============================================================\n');
fprintf('VERIFICA NOMINAL STABILITY (NS)\n');
fprintf('============================================================\n');

poles_nominal = pole(Tmu_scaled);
maxRePole = max(real(poles_nominal));

fprintf('Max Re(polo closed-loop nominale) = %.6e\n', ...
    maxRePole);

if maxRePole < 0

    fprintf('NS = OK: sistema nominalmente stabile.\n');

else

    fprintf('NS = ATTENZIONE: sistema nominalmente instabile.\n');

end

%% ========================================================================
% 12. VERIFICA ROBUST STABILITY (RS)
% ========================================================================

fprintf('\n============================================================\n');
fprintf('VERIFICA ROBUST STABILITY (RS)\n');
fprintf('============================================================\n');

% Modello fisico incerto utilizzato per la sintesi mu
Gmu_phys = G_uncertain_lumped;

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

save('MU_controller_Helicopter.mat', ...
    'K_mu', ...
    'K_mu_scaled', ...
    'CLperf_mu', ...
    'gammaNom_mu', ...
    'stabmarg', ...
    'info_mu', ...
    'info_RS', ...
    'P_mu');

fprintf('\n============================================================\n');
fprintf('SALVATAGGIO COMPLETATO\n');
fprintf('============================================================\n');

fprintf('File "MU_controller_Helicopter.mat" creato correttamente.\n');