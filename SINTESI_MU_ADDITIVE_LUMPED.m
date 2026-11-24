%% ========================================================================
% SINTESI MU SUL MODELLO LUMPED ADDITIVO
%
% Questo script è SEPARATO dalla sintesi moltiplicativa già esistente.
% Legge i risultati precedenti, li usa come riferimento/inizializzazione,
% e salva il nuovo risultato in MU_controller_ADDITIVE.mat.
%
% NON sovrascrive MU_controller.mat né ACTUATOR_LUMPED.mat.
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' MU-SYNTHESIS - MODELLO LUMPED ADDITIVO\n');
fprintf('============================================================\n');

%% 0. CARICAMENTO MODELLI E RISULTATI PRECEDENTI

if ~isfile('HINF_workspace.mat')
    error('HINF_workspace.mat non trovato.');
end
if ~isfile('ACTUATOR_LUMPED_ADDITIVE.mat')
    error('ACTUATOR_LUMPED_ADDITIVE.mat non trovato. Eseguire prima ACTUATOR_LUMPED_ADDITIVE.m.');
end
if ~isfile('MU_controller.mat')
    error('MU_controller.mat non trovato: serve il controllore lumped moltiplicativo già ottenuto.');
end

H = load('HINF_workspace.mat','WS','WU','WT','Du','Dy_inv','G_unc_mu_scaled');
A = load('ACTUATOR_LUMPED_ADDITIVE.mat','G_uncertain_additive_lumped_scaled');
Rold = load('MU_controller.mat','K_mu_scaled_lumped','CLperf_mu_lumped');

if ~isfield(Rold,'K_mu_scaled_lumped') || isempty(Rold.K_mu_scaled_lumped)
    error('K_mu_scaled_lumped non trovato in MU_controller.mat.');
end
K_L = ss(Rold.K_mu_scaled_lumped);

G_add = A.G_uncertain_additive_lumped_scaled;
G_add.InputName = {'ubar1';'ubar2'};
G_add.OutputName = {'alpha_bar';'beta_bar'};

G_NL = H.G_unc_mu_scaled;
G_NL.InputName = {'ubar1';'ubar2'};
G_NL.OutputName = {'alpha_bar';'beta_bar'};

%% 1. PESI E INTERCONNESSIONI GENERALIZZATE

WS_add = H.WS;
WS_add.InputName = {'e1';'e2'};
WS_add.OutputName = {'zS1';'zS2'};

WU_add = H.WU;
WU_add.InputName = {'ubar1';'ubar2'};
WU_add.OutputName = {'zU1';'zU2'};

WT_add = H.WT;
WT_add.InputName = {'alpha_bar';'beta_bar'};
WT_add.OutputName = {'zT1';'zT2'};

SumE1 = sumblk('e1 = r1 - alpha_bar');
SumE2 = sumblk('e2 = r2 - beta_bar');

P_mu_additive = connect( ...
    G_add,WS_add,WU_add,WT_add,SumE1,SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2','e1','e2'});

% Plant generalizzato non-lumped costruito con gli stessi pesi e segnali.
WS_NL = H.WS;
WS_NL.InputName = {'e1';'e2'};
WS_NL.OutputName = {'zS1';'zS2'};
WU_NL = H.WU;
WU_NL.InputName = {'ubar1';'ubar2'};
WU_NL.OutputName = {'zU1';'zU2'};
WT_NL = H.WT;
WT_NL.InputName = {'alpha_bar';'beta_bar'};
WT_NL.OutputName = {'zT1';'zT2'};

P_mu_nonlumped = connect( ...
    G_NL,WS_NL,WU_NL,WT_NL,SumE1,SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2','e1','e2'});

nmeas = 2;
ncont = 2;

%% 2. RIFERIMENTO: VECCHIO CONTROLLore LUMPED MOLTIPLICATIVO

fprintf('\nControllore precedente, sintetizzato sul lumped moltiplicativo:\n');
fprintf('Ordine K_L = %d\n',order(K_L));
if isfield(Rold,'CLperf_mu_lumped')
    fprintf('Performance originaria di sintesi = %.6f\n',Rold.CLperf_mu_lumped);
end

CL_KL_add = lft(P_mu_additive,K_L);
CL_KL_NL = lft(P_mu_nonlumped,K_L);

gamma_KL_add = musynperf(CL_KL_add);
gamma_KL_NL = musynperf(CL_KL_NL);

fprintf('K_L sul nuovo plant additivo: lower = %.6f | upper = %.6f\n', ...
    gamma_KL_add.LowerBound,gamma_KL_add.UpperBound);
fprintf('K_L sul plant non-lumped:    lower = %.6f | upper = %.6f\n', ...
    gamma_KL_NL.LowerBound,gamma_KL_NL.UpperBound);

%% 3. NUOVA SINTESI SUL LUMPED ADDITIVO

% Impostazioni computazionalmente più leggere. Il vecchio controllore viene
% usato solo come inizializzazione; non viene modificato/sovrascritto.
optsMU_additive = musynOptions( ...
    'Display','short', ...
    'MixedMU','on', ...
    'FullDG',true, ...
    'MaxIter',40, ...
    'TolPerf',0.01);

fprintf('\n============================================================\n');
fprintf('AVVIO MUSYN SUL LUMPED ADDITIVO\n');
fprintf('============================================================\n');
fprintf(['Sintesi additiva cold start: il controllore lumped ' ...
    'moltiplicativo è usato solo come riferimento.\n']);
% Sintesi additiva PARTENDO DA ZERO
rng(10);

[K_mu_scaled_additive, CLperf_mu_additive, info_mu_additive] = ...
    musyn(P_mu_additive, nmeas, ncont, optsMU_additive);

% Conserva il controllore originale prodotto da musyn.
K_mu_scaled_additive_full = ss(K_mu_scaled_additive);
N_add_full = order(K_mu_scaled_additive_full);

% Crea il candidato ridotto.
K_mu_scaled_additive_minreal = ...
    minreal(K_mu_scaled_additive_full, 1e-7);

N_add_minreal = order(K_mu_scaled_additive_minreal);

fprintf('\n============================================================\n');
fprintf('MINREAL CONTROLLORE ADDITIVO\n');
fprintf('============================================================\n');

fprintf('Ordine originale   = %d\n', N_add_full);
fprintf('Ordine dopo minreal = %d\n', N_add_minreal);
fprintf('Stati eliminati    = %d\n', N_add_full - N_add_minreal);

fprintf('Performance musyn sul controllore full-order = %.6f\n', ...
    CLperf_mu_additive);

% Da questo punto, le verifiche successive usano il controllore ridotto.
K_mu_scaled_additive = K_mu_scaled_additive_minreal;

fprintf('Ordine del controllore da verificare = %d\n', ...
    order(K_mu_scaled_additive));

%% 4. VALIDAZIONE CON MUSYNPERF

CL_add_synth = lft(P_mu_additive,K_mu_scaled_additive);
CL_add_on_NL = lft(P_mu_nonlumped,K_mu_scaled_additive);

gamma_add_on_add = musynperf(CL_add_synth);
gamma_add_on_NL = musynperf(CL_add_on_NL);

fprintf('\n============================================================\n');
fprintf('RISULTATI DEL CONTROLLore ADDITIVO\n');
fprintf('============================================================\n');
fprintf('K additivo sul lumped additivo: lower = %.6f | upper = %.6f\n', ...
    gamma_add_on_add.LowerBound,gamma_add_on_add.UpperBound);
fprintf('K additivo sul plant non-lumped: lower = %.6f | upper = %.6f\n', ...
    gamma_add_on_NL.LowerBound,gamma_add_on_NL.UpperBound);

if gamma_add_on_add.UpperBound < 1
    fprintf('RP sul nuovo modello lumped additivo: PASS (upper < 1).\n');
else
    fprintf('RP sul nuovo modello lumped additivo: NON certificata (upper >= 1).\n');
end
if gamma_add_on_NL.UpperBound < 1
    fprintf('RP sul modello non-lumped: PASS (upper < 1).\n');
else
    fprintf('RP sul modello non-lumped: NON certificata (upper >= 1).\n');
end

%% ========================================================================
% 4.1 WCGain: CONFRONTO SULLO STESSO PLANT
% ========================================================================

fprintf('\n============================================================\n');
fprintf('VALIDAZIONE ROBUST PERFORMANCE CON WCGain\n');
fprintf('============================================================\n');

% Stessa opzione di analisi per tutti i casi.
% 'a' forza la stima del limite superiore di mu mediante LMI.
optsWC = wcOptions('MussvOptions','a');

% ------------------------------------------------------------
% Controllore di riferimento: lumped moltiplicativo
% ------------------------------------------------------------

fprintf('\nWCGain: K_L sul plant additivo...\n');
[wc_KL_add, wcu_KL_add] = wcgain(CL_KL_add,optsWC);

fprintf('WCGain: K_L sul plant non-lumped...\n');
[wc_KL_NL, wcu_KL_NL] = wcgain(CL_KL_NL,optsWC);

% ------------------------------------------------------------
% Nuovo controllore: lumped additivo
% ------------------------------------------------------------

fprintf('\nWCGain: K additivo sul plant additivo...\n');
[wc_add_on_add, wcu_add_on_add] = ...
    wcgain(CL_add_synth,optsWC);

fprintf('WCGain: K additivo sul plant non-lumped...\n');
[wc_add_on_NL, wcu_add_on_NL] = ...
    wcgain(CL_add_on_NL,optsWC);

fprintf('\n------------------------------------------------------------\n');
fprintf('RISULTATI WCGain - LOWER / UPPER\n');
fprintf('------------------------------------------------------------\n');

fprintf('K_L sul plant additivo:   %.6f / %.6f\n', ...
    wc_KL_add.LowerBound,wc_KL_add.UpperBound);

fprintf('K_L sul plant non-lumped: %.6f / %.6f\n', ...
    wc_KL_NL.LowerBound,wc_KL_NL.UpperBound);

fprintf('K_A sul plant additivo:   %.6f / %.6f\n', ...
    wc_add_on_add.LowerBound,wc_add_on_add.UpperBound);

fprintf('K_A sul plant non-lumped: %.6f / %.6f\n', ...
    wc_add_on_NL.LowerBound,wc_add_on_NL.UpperBound);

% Esito RP: conta il limite superiore.
if wc_add_on_add.UpperBound < 1
    fprintf('RP sul plant additivo: PASS\n');
else
    fprintf('RP sul plant additivo: NON CERTIFICATA\n');
end

if wc_add_on_NL.UpperBound < 1
    fprintf('RP sul plant non-lumped: PASS\n');
else
    fprintf('RP sul plant non-lumped: NON CERTIFICATA\n');
end


%% ========================================================================
% 4.2 NOMINAL STABILITY E ROBUST STABILITY SUL NON-LUMPED
% ========================================================================

fprintf('\n============================================================\n');
fprintf('NOMINAL STABILITY E ROBUST STABILITY\n');
fprintf('============================================================\n');

% Plant non-lumped normalizzato, senza pesi di prestazione.
G_NL_RS = H.G_unc_mu_scaled;
G_NL_RS.InputName  = {'u1';'u2'};
G_NL_RS.OutputName = {'y1';'y2'};

% Connessioni di feedback: e = r - y.
SumRS1 = sumblk('e1 = r1 - y1');
SumRS2 = sumblk('e2 = r2 - y2');

% Closed-loop del controllore additivo.
K_RS_A = ss(K_mu_scaled_additive);
K_RS_A.InputName  = {'e1';'e2'};
K_RS_A.OutputName = {'u1';'u2'};

CL_RS_add_NL = connect( ...
    G_NL_RS,K_RS_A,SumRS1,SumRS2, ...
    {'r1','r2'},{'y1','y2'});

% Verifica nominale.
CL_RS_add_NL_nom = getNominal(CL_RS_add_NL);

nominalStableAdd = isstable(CL_RS_add_NL_nom);

fprintf('Stabilita'' nominale K_A sul non-lumped = %d\n', ...
    nominalStableAdd);

% Verifica robusta.
[stabmarg_add_NL,destabunc_add_NL,info_RS_add_NL] = ...
    robuststab(CL_RS_add_NL);

fprintf('RS lower bound K_A sul non-lumped = %.6f\n', ...
    stabmarg_add_NL.LowerBound);

fprintf('RS upper bound K_A sul non-lumped = %.6f\n', ...
    stabmarg_add_NL.UpperBound);

if stabmarg_add_NL.LowerBound > 1
    fprintf('RS sul non-lumped: PASS (lower bound > 1).\n');
else
    fprintf('RS sul non-lumped: NON CERTIFICATA.\n');
end

%% ========================================================================
% 4.3 CONFRONTO FULL-ORDER vs MINREAL
% ========================================================================

fprintf('\n============================================================\n');
fprintf('CONFRONTO CONTROLLORE FULL-ORDER E MINREAL\n');
fprintf('============================================================\n');

% Il controllore attualmente verificato è quello MINREAL.
% Qui si verifica separatamente il full-order originario.

K_full_test = ss(K_mu_scaled_additive_full);

fprintf('Ordine full-order = %d\n',order(K_full_test));
fprintf('Ordine minreal    = %d\n',order(K_mu_scaled_additive_minreal));

fprintf('Differenza ingresso-uscita full/minreal:\n');
fprintf('||K_full-K_minreal||inf = %.6e\n', ...
    norm(K_mu_scaled_additive_full - ...
    K_mu_scaled_additive_minreal,inf));

% Closed-loop pesati sul lumped additivo
CL_full_add = lft(P_mu_additive,K_full_test);

% Closed-loop pesato sul non-lumped
CL_full_NL = lft(P_mu_nonlumped,K_full_test);

% Verifica RP con musynperf
gamma_full_add = musynperf(CL_full_add);
gamma_full_NL = musynperf(CL_full_NL);

fprintf('\nMUSYNPERF FULL-ORDER:\n');
fprintf('Lumped additivo: lower = %.6f | upper = %.6f\n', ...
    gamma_full_add.LowerBound,gamma_full_add.UpperBound);

fprintf('Non-lumped:      lower = %.6f | upper = %.6f\n', ...
    gamma_full_NL.LowerBound,gamma_full_NL.UpperBound);

% Verifica RP con wcgain
fprintf('\nWCGain FULL-ORDER:\n');

[wc_full_add,wcu_full_add] = wcgain(CL_full_add,optsWC);
[wc_full_NL,wcu_full_NL] = wcgain(CL_full_NL,optsWC);

fprintf('Lumped additivo: lower = %.6f | upper = %.6f\n', ...
    wc_full_add.LowerBound,wc_full_add.UpperBound);

fprintf('Non-lumped:      lower = %.6f | upper = %.6f\n', ...
    wc_full_NL.LowerBound,wc_full_NL.UpperBound);

% Verifica RS sul modello non-lumped
K_RS_full = K_full_test;
K_RS_full.InputName  = {'e1';'e2'};
K_RS_full.OutputName = {'u1';'u2'};

CL_RS_full_NL = connect( ...
    G_NL_RS,K_RS_full,SumRS1,SumRS2, ...
    {'r1','r2'},{'y1','y2'});

nominalStableFull = isstable(getNominal(CL_RS_full_NL));

[stabmarg_full_NL,destabunc_full_NL,info_RS_full_NL] = ...
    robuststab(CL_RS_full_NL);

fprintf('\nROBUST STABILITY FULL-ORDER:\n');
fprintf('Stabilita'' nominale = %d\n',nominalStableFull);
fprintf('RS lower = %.6f | upper = %.6f\n', ...
    stabmarg_full_NL.LowerBound, ...
    stabmarg_full_NL.UpperBound);

if wc_full_NL.UpperBound < 1 && ...
        stabmarg_full_NL.LowerBound > 1 && nominalStableFull
    fprintf('FULL-ORDER: RP e RS certificate sul non-lumped.\n');
else
    fprintf('FULL-ORDER: almeno una verifica sul non-lumped non è certificata.\n');
end

fprintf('\nRISULTATI MINREAL GIA'' CALCOLATI:\n');
fprintf('RP upper wcgain non-lumped = %.6f\n', ...
    wc_add_on_NL.UpperBound);
fprintf('RS lower non-lumped       = %.6f\n', ...
    stabmarg_add_NL.LowerBound);

%% ========================================================================
% 4.4 SELEZIONE DEL CONTROLLORE VALIDATO
% ========================================================================

pass_minreal = ...
    nominalStableAdd && ...
    wc_add_on_add.UpperBound < 1 && ...
    wc_add_on_NL.UpperBound < 1 && ...
    stabmarg_add_NL.LowerBound > 1;

pass_full = ...
    nominalStableFull && ...
    wc_full_add.UpperBound < 1 && ...
    wc_full_NL.UpperBound < 1 && ...
    stabmarg_full_NL.LowerBound > 1;

if pass_minreal

    K_mu_scaled_additive = K_mu_scaled_additive_minreal;
    controllerSelected = 'MINREAL';

    fprintf('\nSelezionato MINREAL: RP e RS certificate.\n');

elseif pass_full

    K_mu_scaled_additive = K_mu_scaled_additive_full;
    controllerSelected = 'FULL_ORDER';

    fprintf('\nSelezionato FULL-ORDER: RP e RS certificate.\n');

else

    controllerSelected = 'NONE_CERTIFIED';
    K_mu_scaled_additive = [];

    warning(['Nessuno dei due controllori ha superato tutte le ' ...
        'verifiche. Nessun controllore finale selezionato.']);

end

fprintf('Controllore selezionato = %s\n',controllerSelected);

%% 5. DESCALATURA E SALVATAGGIO SEPARATO

if strcmp(controllerSelected,'NONE_CERTIFIED')
    K_mu_additive_phys = [];
else
    K_mu_additive_phys = ...
        H.Du*K_mu_scaled_additive*H.Dy_inv;

    K_mu_additive_phys.InputName = {'e1','e2'};
    K_mu_additive_phys.OutputName = {'u1','u2'};
end

outFile = 'MU_controller_ADDITIVE_15PCT_COLDSTART_MINREAL.mat';
if isfile(outFile)
    backupName = ['MU_controller_ADDITIVE_previous_' ...
        char(datetime('now','Format','yyyyMMdd_HHmmss')) '.mat'];
    copyfile(outFile,backupName);
    fprintf('Versione additiva precedente salvata in %s\n',backupName);
end

save(outFile, ...
    'K_mu_scaled_additive','K_mu_additive_phys', ...
    'K_mu_scaled_additive_full', ...
    'K_mu_scaled_additive_minreal', ...
    'N_add_full', ...
    'N_add_minreal', ...
    'CLperf_mu_additive','info_mu_additive','optsMU_additive', ...
    'gamma_KL_add','gamma_KL_NL', ...
    'gamma_add_on_add','gamma_add_on_NL', ...
    'P_mu_additive','P_mu_nonlumped', ...
    'optsWC', ...
    'gamma_full_add','gamma_full_NL', ...
    'wc_full_add','wcu_full_add', ...
    'wc_full_NL','wcu_full_NL', ...
    'nominalStableFull', ...
    'CL_RS_full_NL', ...
    'stabmarg_full_NL','destabunc_full_NL','info_RS_full_NL', ...
    'wc_KL_add','wcu_KL_add', ...
    'wc_KL_NL','wcu_KL_NL', ...
    'wc_add_on_add','wcu_add_on_add', ...
    'wc_add_on_NL','wcu_add_on_NL', ...
    'CL_RS_add_NL','CL_RS_add_NL_nom', ...
    'nominalStableAdd', ...
    'controllerSelected', ...
    'stabmarg_add_NL','destabunc_add_NL','info_RS_add_NL');

fprintf('\n============================================================\n');
fprintf('NUOVI RISULTATI SALVATI IN %s\n',outFile);
fprintf('MU_controller.mat non è stato modificato.\n');
fprintf('============================================================\n');
