%% ========================================================================
% MU-SYNTHESIS - LUMPED MOLTIPLICATIVO CON TAU_D ESPLICITO
% ESECUZIONE NUOVA DA CAPO + VERIFICHE NS, NP, RS, RP
%
% - Ricostruisce il plant generalizzato usando il modello tau esplicito.
% - Riparte dal controllore lumped moltiplicativo originale (Kinit), quindi
%   la D-K iteration ricomincia dall'iterazione 1.
% - NON usa info_tau_explicit(end).K e NON continua la run precedente.
% - NON sovrascrive MU_controller.mat né i risultati tau esplicito precedenti.
% - Salva prima il candidato grezzo e poi il risultato con le verifiche.
% ========================================================================

close all;
clc;

fprintf('============================================================\n');
fprintf(' MU-SYNTHESIS - TAU_D ESPLICITO\n');
fprintf('============================================================\n');

%% 1. CONTROLLI E CARICAMENTO DEI FILE

assert(isfile('HINF_workspace.mat'), ...
    'HINF_workspace.mat non trovato nella cartella corrente.');
assert(isfile('ACTUATOR_LUMPED_TAU_ESPLICITO.mat'), ...
    'Eseguire prima ACTUATOR_LUMPED_TAU_ESPLICITO.m.');
assert(isfile('MU_controller.mat'), ...
    'MU_controller.mat non trovato: serve il K lumped moltiplicativo iniziale.');

H = load('HINF_workspace.mat','WS','WU','WT');
A = load('ACTUATOR_LUMPED_TAU_ESPLICITO.mat', ...
    'G_uncertain_lumped_tau_explicit_scaled');
Kdata = load('MU_controller.mat', ...
    'K_mu_scaled_lumped','CLperf_mu_lumped');

assert(isfield(Kdata,'K_mu_scaled_lumped') && ...
       ~isempty(Kdata.K_mu_scaled_lumped), ...
    'K_mu_scaled_lumped non disponibile in MU_controller.mat.');

G = A.G_uncertain_lumped_tau_explicit_scaled;
Kinit = ss(Kdata.K_mu_scaled_lumped);

fprintf('Ordine K iniziale = %d\n',order(Kinit));
fprintf('Performance lumped precedente = %.6f\n', ...
    Kdata.CLperf_mu_lumped);

%% 2. COSTRUZIONE DEL PLANT GENERALIZZATO

G.InputName = {'ubar1','ubar2'};
G.OutputName = {'alpha_bar','beta_bar'};

WS = ss(H.WS);
WS.InputName = {'e1','e2'};
WS.OutputName = {'zS1','zS2'};

WU = ss(H.WU);
WU.InputName = {'ubar1','ubar2'};
WU.OutputName = {'zU1','zU2'};

WT = ss(H.WT);
WT.InputName = {'alpha_bar','beta_bar'};
WT.OutputName = {'zT1','zT2'};

SumE1 = sumblk('e1 = r1 - alpha_bar');
SumE2 = sumblk('e2 = r2 - beta_bar');

P_mu_tau_explicit = connect( ...
    G,WS,WU,WT,SumE1,SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2','e1','e2'});

fprintf('\nPlant generalizzato tau esplicito creato.\n');

%% 3. VERIFICA PRE-SINTESI DEL CONTROLLore INIZIALE

CLinit_nom = lft(P_mu_tau_explicit.NominalValue,Kinit);
CLinit_unc = lft(P_mu_tau_explicit,Kinit);

NS_init = isstable(CLinit_nom);
fprintf('Stabilita nominale K iniziale = %d\n',NS_init);
if ~NS_init
    error(['Il controllore lumped iniziale non stabilizza nominalmente ' ...
        'il plant tau esplicito. Sintesi interrotta.']);
end

perf_init = musynperf(CLinit_unc);
fprintf('\n============================================\n');
fprintf('PRESTAZIONE INIZIALE K LUMPED\n');
fprintf('============================================\n');
fprintf('musynperf lower bound = %.6f\n',perf_init.LowerBound);
fprintf('musynperf upper bound = %.6f\n',perf_init.UpperBound);

%% 4. OPZIONI D-K PER UNA NUOVA RUN

opts_tau = musynOptions( ...
    'Display','short', ...
    'MixedMU','on', ...
    'FullDG',true, ...
    'MaxIter',40, ...
    'TolPerf',0.01);

fprintf('\n============================================\n');
fprintf('AVVIO MUSYN\n');
fprintf('============================================\n');

rng(10);

[K_mu_scaled_tau_explicit, ...
    CLperf_tau_explicit, ...
    info_tau_explicit] = musyn( ...
    P_mu_tau_explicit,2,2,opts_tau);

K_mu_scaled_tau_explicit = ss(K_mu_scaled_tau_explicit);

fprintf('\n============================================\n');
fprintf('RISULTATO DELLA NUOVA SINTESI TAU ESPLICITO\n');
fprintf('============================================\n');
fprintf('CLperf musyn = %.6f\n',CLperf_tau_explicit);
fprintf('Ordine K restituito = %d\n',order(K_mu_scaled_tau_explicit));

%% 5. SALVATAGGIO GREZZO IMMEDIATO DOPO MUSYN
% Protegge il controllore e la cronologia anche se un'analisi successiva
% dovesse generare un errore numerico.

save('MU_controller_TAU_ESPLICITO_RERUN_RAW.mat', ...
    'K_mu_scaled_tau_explicit','CLperf_tau_explicit', ...
    'info_tau_explicit','Kinit','perf_init','opts_tau', ...
    'P_mu_tau_explicit','Kdata');

fprintf('Candidato grezzo salvato in MU_controller_TAU_ESPLICITO_RERUN_RAW.mat\n');

%% 6. VERIFICHE NS, NP, RS, RP PER K INIZIALE E CANDIDATO

Klist = {Kinit, K_mu_scaled_tau_explicit};
labels = {'K iniziale lumped', 'K candidato nuova sintesi'};

NS_all = false(1,2);
NP_all = Inf(1,2);
fPeakNP_all = NaN(1,2);
RS_lower_all = NaN(1,2);
RS_upper_all = NaN(1,2);
RP_lower_all = NaN(1,2);
RP_upper_all = NaN(1,2);

CL_unc_all = cell(1,2);
CL_nom_all = cell(1,2);
stabmarg_all = cell(1,2);
destabunc_all = cell(1,2);
reportRS_all = cell(1,2);
infoRS_all = cell(1,2);
wc_all = cell(1,2);
wcu_all = cell(1,2);
perf_all = cell(1,2);

optsRS_tau = robuststabOptions('Display','off');
optsWC_tau = wcOptions('Display','off','MussvOptions','a');

for i = 1:2
    fprintf('\n============================================================\n');
    fprintf('VERIFICHE: %s\n',labels{i});
    fprintf('Ordine K = %d\n',order(Klist{i}));
    fprintf('============================================================\n');

    CL_nom_all{i} = lft(P_mu_tau_explicit.NominalValue,Klist{i});
    CL_unc_all{i} = lft(P_mu_tau_explicit,Klist{i});

    %% NS - Nominal Stability
    NS_all(i) = isstable(CL_nom_all{i});
    fprintf('\nNS (Nominal Stability) = %d\n',NS_all(i));
    if NS_all(i)
        fprintf('NS PASS: closed-loop nominalmente stabile.\n');
    else
        fprintf('NS FAIL: closed-loop nominalmente instabile.\n');
    end

    %% NP - Nominal Performance
    if NS_all(i)
        [NP_all(i),fPeakNP_all(i)] = hinfnorm(CL_nom_all{i});
    else
        NP_all(i) = Inf;
        fPeakNP_all(i) = NaN;
    end
    fprintf('\nNP gamma = %.6f\n',NP_all(i));
    fprintf('Frequenza di picco NP = %.6f rad/s\n',fPeakNP_all(i));
    if NP_all(i) < 1
        fprintf('NP PASS: gamma nominale < 1.\n');
    else
        fprintf('NP FAIL: gamma nominale >= 1 oppure instabilita nominale.\n');
    end

    %% RS - Robust Stability
    [stabmarg_all{i},destabunc_all{i},reportRS_all{i},infoRS_all{i}] = ...
        robuststab(CL_unc_all{i},optsRS_tau);

    RS_lower_all(i) = stabmarg_all{i}.LowerBound;
    RS_upper_all(i) = stabmarg_all{i}.UpperBound;
    fprintf('\nRS lower bound = %.6f\n',RS_lower_all(i));
    fprintf('RS upper bound = %.6f\n',RS_upper_all(i));
    if RS_lower_all(i) > 1
        fprintf('RS PASS: stabilita robusta certificata.\n');
    elseif RS_upper_all(i) < 1
        fprintf('RS FAIL: il margine indica instabilita entro l''incertezza specificata.\n');
    else
        fprintf('RS INCONCLUSIVA: i bound non certificano un esito.\n');
    end

    %% RP - Robust Performance mediante wcgain
    [wc_all{i},wcu_all{i}] = wcgain(CL_unc_all{i},optsWC_tau);
    RP_lower_all(i) = wc_all{i}.LowerBound;
    RP_upper_all(i) = wc_all{i}.UpperBound;
    fprintf('\nRP wcgain lower bound = %.6f\n',RP_lower_all(i));
    fprintf('RP wcgain upper bound = %.6f\n',RP_upper_all(i));
    if isfield(wc_all{i},'CriticalFrequency')
        fprintf('Frequenza critica RP = %.6f rad/s\n', ...
            wc_all{i}.CriticalFrequency);
    end
    if RP_upper_all(i) < 1
        fprintf('RP PASS: prestazione robusta certificata (upper bound < 1).\n');
    elseif RP_lower_all(i) > 1
        fprintf('RP FAIL: prestazione robusta sicuramente > 1.\n');
    else
        fprintf('RP INCONCLUSIVA: i bound comprendono 1.\n');
    end

    %% Analisi supplementare musynperf
    perf_all{i} = musynperf(CL_unc_all{i});
    fprintf('\nmusynperf lower bound = %.6f\n',perf_all{i}.LowerBound);
    fprintf('musynperf upper bound = %.6f\n',perf_all{i}.UpperBound);
end

%% 7. ALIAS DEI RISULTATI DEL CANDIDATO PER ACCESSO SEMPLICE

NS_tau = NS_all(2);
NP_tau = NP_all(2);
fPeakNP_tau = fPeakNP_all(2);
stabmarg_tau = stabmarg_all{2};
destabunc_tau = destabunc_all{2};
reportRS_tau = reportRS_all{2};
infoRS_tau = infoRS_all{2};
RS_lower_tau = RS_lower_all(2);
RS_upper_tau = RS_upper_all(2);
wc_tau = wc_all{2};
wcu_tau = wcu_all{2};
RP_lower_tau = RP_lower_all(2);
RP_upper_tau = RP_upper_all(2);
CL_tau_unc = CL_unc_all{2};
CL_tau_nom = CL_nom_all{2};
perf_final_tau = perf_all{2};

%% 8. CONFRONTO RIASSUNTIVO

fprintf('\n============================================================\n');
fprintf('CONFRONTO FINALE\n');
fprintf('============================================================\n');
fprintf('%-28s %12s %12s\n','Metrica','K iniziale','K candidato');
fprintf('%-28s %12d %12d\n','Ordine',order(Klist{1}),order(Klist{2}));
fprintf('%-28s %12d %12d\n','NS',NS_all(1),NS_all(2));
fprintf('%-28s %12.6f %12.6f\n','NP gamma',NP_all(1),NP_all(2));
fprintf('%-28s %12.6f %12.6f\n','RS lower',RS_lower_all(1),RS_lower_all(2));
fprintf('%-28s %12.6f %12.6f\n','RS upper',RS_upper_all(1),RS_upper_all(2));
fprintf('%-28s %12.6f %12.6f\n','RP wcgain lower',RP_lower_all(1),RP_lower_all(2));
fprintf('%-28s %12.6f %12.6f\n','RP wcgain upper',RP_upper_all(1),RP_upper_all(2));
fprintf('%-28s %12.6f %12.6f\n','musynperf upper', ...
    perf_all{1}.UpperBound,perf_all{2}.UpperBound);

%% 9. SALVATAGGIO COMPLETO SEPARATO
% Non sovrascrive MU_controller.mat né MU_controller_TAU_ESPLICITO.mat.

save('MU_controller_TAU_ESPLICITO_RERUN.mat', ...
    'K_mu_scaled_tau_explicit','CLperf_tau_explicit','info_tau_explicit', ...
    'perf_init','perf_final_tau','Kinit','opts_tau','P_mu_tau_explicit', ...
    'NS_tau','NP_tau','fPeakNP_tau','stabmarg_tau','destabunc_tau', ...
    'reportRS_tau','infoRS_tau','RS_lower_tau','RS_upper_tau', ...
    'wc_tau','wcu_tau','RP_lower_tau','RP_upper_tau', ...
    'CL_tau_unc','CL_tau_nom','NS_all','NP_all','fPeakNP_all', ...
    'RS_lower_all','RS_upper_all','RP_lower_all','RP_upper_all', ...
    'CL_unc_all','CL_nom_all','stabmarg_all','destabunc_all', ...
    'reportRS_all','infoRS_all','wc_all','wcu_all','perf_all','labels', ...
    'Kdata');

fprintf('\n============================================================\n');
fprintf('SALVATAGGIO COMPLETATO\n');
fprintf('Candidato grezzo: MU_controller_TAU_ESPLICITO_RERUN_RAW.mat\n');
fprintf('Risultati e verifiche: MU_controller_TAU_ESPLICITO_RERUN.mat\n');
fprintf('I file di sintesi precedenti non sono stati modificati.\n');
