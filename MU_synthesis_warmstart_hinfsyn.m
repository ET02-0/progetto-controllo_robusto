%% ========================================================================
% SINTESI H-INFINITO MU-SYNTHESIS (D-K ITERATION) - ELICOTTERO 2-DOF
% ========================================================================

close all;
clc;

%% INTERRUTTORI DI ESECUZIONE
% true: esegue anche la sintesi diretta sul modello non-lumped.
% false: salta la seconda sintesi ma conserva la validazione del K lumped.
RUN_MUSYN_NONLUMPED = true;

% true: cerca un controllore di ordine ridotto e lo verifica sul non-lumped.
% false: mantiene il controllore full-order senza effettuare la ricerca.
RUN_ORDER_REDUCTION = false;

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
%% ========================================================================
% CARICAMENTO DEI RISULTATI DELLA PRECEDENTE SINTESI MU
% ========================================================================

resumeFile = 'MU_controller.mat';

if ~isfile(resumeFile)
    error(['File %s non trovato. Serve il MAT della precedente ' ...
        'esecuzione, contenente info_mu_lumped e info_mu_nonlumped.'], ...
        resumeFile);
end

resumeData = load(resumeFile);

% Indici utilizzati SOLO alla prima ripresa

resumeIterNonlumped = 8;

%% LUMPED: CARICAMENTO DEL CONTROLLore GIA' VALIDATO
% Non viene eseguita una nuova sintesi lumped.
% Si riutilizza il miglior controllore lumped precedentemente salvato.

if isfield(resumeData,'K_mu_scaled_lumped_resume') && ...
        isfield(resumeData,'CLperf_mu_lumped_resume') && ...
        ~isempty(resumeData.K_mu_scaled_lumped_resume)

    K_mu_scaled_lumped = ss(resumeData.K_mu_scaled_lumped_resume);
    CLperf_mu_lumped = resumeData.CLperf_mu_lumped_resume;

elseif isfield(resumeData,'K_mu_scaled_lumped') && ...
        isfield(resumeData,'CLperf_mu_lumped')

    % Compatibilità con salvataggi precedenti senza variabili _resume.
    K_mu_scaled_lumped = ss(resumeData.K_mu_scaled_lumped);
    CLperf_mu_lumped = resumeData.CLperf_mu_lumped;

else
    error(['Non trovo un controllore lumped salvato. ' ...
        'Impossibile procedere senza una precedente sintesi lumped.']);
end

if isfield(resumeData,'info_mu_lumped')
    info_mu_lumped = resumeData.info_mu_lumped;
else
    info_mu_lumped = [];
end

fprintf('\nLUMPED: nuova sintesi disattivata.\n');
fprintf('Riutilizzo il controllore lumped precedentemente salvato.\n');
fprintf('Ordine K lumped = %d\n',order(K_mu_scaled_lumped));
fprintf('Performance di sintesi lumped = %.6f\n',CLperf_mu_lumped);

%% NON-LUMPED: CARICAMENTO DEL PRECEDENTE RISULTATO
% Il controllore precedente viene mantenuto come benchmark.


if RUN_MUSYN_NONLUMPED

    if isfield(resumeData,'K_mu_scaled_nonlumped_resume') && ...
            isfield(resumeData,'CLperf_mu_nonlumped_resume') && ...
            ~isempty(resumeData.K_mu_scaled_nonlumped_resume)

        Kprev_mu_nonlumped = ...
            ss(resumeData.K_mu_scaled_nonlumped_resume);

        CLperf_prev_nonlumped = ...
            resumeData.CLperf_mu_nonlumped_resume;

    elseif isfield(resumeData,'K_mu_scaled_nonlumped') && ...
            isfield(resumeData,'CLperf_mu_nonlumped')

        % Compatibilità con salvataggi precedenti senza variabili _resume.
        Kprev_mu_nonlumped = ...
            ss(resumeData.K_mu_scaled_nonlumped);

        CLperf_prev_nonlumped = ...
            resumeData.CLperf_mu_nonlumped;

    else
        error(['Non trovo il precedente controllore non-lumped ' ...
            'da usare come benchmark.']);
    end

    fprintf('\nPrecedente controllore non-lumped caricato.\n');
    fprintf('Ordine = %d\n',order(Kprev_mu_nonlumped));
    fprintf('Performance = %.6f\n',CLperf_prev_nonlumped);

    if isfield(resumeData,'info_mu_nonlumped')
        info_mu_nonlumped = resumeData.info_mu_nonlumped;
    else
        info_mu_nonlumped = [];
    end

end

fprintf('\nConfigurazione della sintesi:\n');
fprintf('  LUMPED: controllore esistente, nessuna nuova sintesi\n');

if RUN_MUSYN_NONLUMPED
    fprintf('  NON-LUMPED: nuova sintesi abilitata\n');
end

requiredNonLumped = {
    'J_alpha'
    'l'
    'omega_n1'
    'omega_n2'
    'tau_d1'
    'tau_d2'
    'Delta_Pade1'
    'Delta_Pade2'
};

availableUncertainty = ...
    fieldnames(G_unc_mu_scaled.Uncertainty);

missingUncertainty = ...
    setdiff(requiredNonLumped, availableUncertainty);

if ~isempty(missingUncertainty)
    error(['G_unc_mu_scaled non contiene tutte le incertezze ' ...
           'attese. Rieseguire HINF_SETUP. Mancano: %s'], ...
           strjoin(missingUncertainty, ', '));
end

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
% 5.1 VERIFICA E IMPIANTO GENERALIZZATO NON-LUMPED
% ========================================================================

% Controlla che il setup abbia salvato il modello per la sintesi.
if ~exist('G_unc_mu_scaled_synth','var')
    error(['G_unc_mu_scaled_synth non trovato. Rieseguire HINF_SETUP ' ...
        'e verificare che il modello sia incluso in HINF_workspace.mat.']);
end

% Mantieni separati i due modelli:
% - originale: analisi e validazione;
% - semplificato: sintesi MUSYN.

Gmu_nonlumped_validation = G_unc_mu_scaled;
Gmu_nonlumped_synth = G_unc_mu_scaled_synth;

%% VERIFICA DEL MODELLO COMPLETO NON-LUMPED SEMPLIFICATO

Gnom_original = Gmu_nonlumped_validation.NominalValue;
Gnom_synth    = Gmu_nonlumped_synth.NominalValue;

Dnom_final = minreal( ...
    Gnom_original - Gnom_synth, 1e-8);

errNom_final = norm(Dnom_final,inf);

fprintf('\nDifferenza nominale modello completo = %.6e\n', ...
    errNom_final);

% Confronto su campioni comuni delle incertezze.
Ns_final = 50;

[G_original_samples,SampleValues_final] = ...
    usample(Gmu_nonlumped_validation,Ns_final);

G_synth_samples = usubs( ...
    Gmu_nonlumped_synth,SampleValues_final);

err_final = zeros(Ns_final,1);

for k = 1:Ns_final
    err_final(k) = norm( ...
        G_original_samples(:,:,k) - ...
        G_synth_samples(:,:,k),inf);
end

fprintf('Errore massimo sui %d campioni del modello completo = %.6e\n', ...
    Ns_final,max(err_final));

if errNom_final > 1e-6 || max(err_final) > 1e-6
    error(['Il modello completo semplificato non ha superato ' ...
        'il controllo di coerenza. Non procedere con musyn.']);
end

fprintf('Coerenza numerica del modello completo verificata.\n');

%% COSTRUZIONE DEL PLANT GENERALIZZATO PER LA SINTESI

Gmu_nonlumped_synth.InputName = {'ubar1','ubar2'};
Gmu_nonlumped_synth.OutputName = {'alpha_bar','beta_bar'};

P_mu_nonlumped = connect( ...
    Gmu_nonlumped_synth, ...
    WS, ...
    WU, ...
    WT, ...
    SumE1, ...
    SumE2, ...
    {'r1','r2','ubar1','ubar2'}, ...
    {'zS1','zS2', ...
    'zU1','zU2', ...
    'zT1','zT2', ...
    'e1','e2'});

fprintf('\nPlant generalizzato non-lumped creato.\n');

%% ========================================================================
% 6. D-K ITERATION - MUSYN
% ========================================================================

nmeas = 2;
ncont = 2;

% Nuova sintesi NON-LUMPED, inizializzata con il controllore
% LUMPED precedentemente validato sul modello non-lumped.
optsMU_nonlumped = musynOptions( ...
    'Display','short', ...
    'MixedMU','on', ...
    'FullDG',true, ...
    'TargetPerf',0, ...
    'MaxIter',40, ...
    'TolPerf',0.01);

%% RIUTILIZZO DEL CONTROLLore LUMPED ESISTENTE

fprintf('\n============================================================\n');
fprintf('SINTESI LUMPED DISATTIVATA\n');
fprintf('============================================================\n');

% Si mantiene il controllore lumped caricato da MU_controller.mat.
% Non viene chiamato musyn sul plant lumped.

K_mu_scaled_lumped = ss(K_mu_scaled_lumped);

K_mu_scaled = K_mu_scaled_lumped;
CLperf_mu = CLperf_mu_lumped;
info_mu = info_mu_lumped;

fprintf('Controllore lumped mantenuto: ordine %d\n', ...
    order(K_mu_scaled_lumped));

fprintf('Performance lumped salvata = %.6f\n', ...
    CLperf_mu_lumped);
%% VERIFICA DEL CONTROLLore LUMPED SELEZIONATO

Kraw = ss(K_mu_scaled_lumped);

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
% 6.1 D-K ITERATION DIRETTA SUL MODELLO NON-LUMPED (OPZIONALE)
% ========================================================================

if RUN_MUSYN_NONLUMPED

    %% MUSYN DIRETTA NON-LUMPED - KINIT DA MIXSYN

    Hc = load('HINF_controllers.mat');

    if ~isfield(Hc,'K_mix_scaled')
        error('K_mix_scaled non trovato in HINF_controllers.mat.');
    end

    Kinit_mu_nonlumped = ss(Hc.K_mix_scaled);

    % Controllo dimensioni: 2 ingressi di misura, 2 uscite di controllo
    if size(Kinit_mu_nonlumped,1) ~= 2 || ...
            size(Kinit_mu_nonlumped,2) ~= 2
        error('K_mix_scaled non ha dimensioni 2x2.');
    end

    % Verifica che il controllore stabilizzi il plant nominale
    CL_init_nom = lft( ...
        P_mu_nonlumped.NominalValue, ...
        Kinit_mu_nonlumped);

    fprintf('Ordine Kinit mixsyn = %d\n', ...
        order(Kinit_mu_nonlumped));

    fprintf('Stabilita'' closed-loop nominale Kinit = %d\n', ...
        isstable(CL_init_nom));

    if ~isstable(CL_init_nom)
        error(['Il controllore mixsyn non stabilizza il plant nominale. ' ...
            'Non procedere con questa inizializzazione.']);
    end

    %% OPZIONI MUSYN
    optsMU_nonlumped = musynOptions( ...
        'Display','short', ...
        'MixedMU','on', ...
        'FullDG',true, ...
        'FitOrder',[5 2], ...
        'TargetPerf',0, ...
        'MaxIter',40, ...
        'TolPerf',0.01);

    rng(10);

    %% SINTESI
    [K_direct_mixsyn, CLperf_direct_mixsyn, ...
        info_direct_mixsyn] = ...
        musyn(P_mu_nonlumped, nmeas, ncont, ...
        Kinit_mu_nonlumped, optsMU_nonlumped);

    fprintf('\nPerformance MUSYN = %.6f\n', ...
        CLperf_direct_mixsyn);

    fprintf('Ordine K restituito = %d\n', ...
        order(K_direct_mixsyn));

    %% SALVATAGGIO SEPARATO
    save('MU_DIRECT_FROM_MIXSYN_INIT.mat', ...
        'Kinit_mu_nonlumped', ...
        'K_direct_mixsyn', ...
        'CLperf_direct_mixsyn', ...
        'info_direct_mixsyn', ...
        'optsMU_nonlumped');

    %% CONFRONTO CON IL PRECEDENTE CONTROLLore NON-LUMPED

    % La nuova sintesi parte da K_mix_scaled.
    % Il precedente controllore non-lumped viene utilizzato
    % esclusivamente come benchmark per la selezione.

    if CLperf_direct_mixsyn <= CLperf_prev_nonlumped

        % Il nuovo controllore ha una performance migliore
        K_mu_scaled_nonlumped = ss(K_direct_mixsyn);
        CLperf_mu_nonlumped = CLperf_direct_mixsyn;
        info_mu_nonlumped = info_direct_mixsyn;

        fprintf(['Nuovo controllore non-lumped selezionato: ' ...
            'performance migliorata o confermata.\n']);

    else

        % Il nuovo controllore è peggiore:
        % si mantiene quello precedentemente salvato.
        K_mu_scaled_nonlumped = Kprev_mu_nonlumped;
        CLperf_mu_nonlumped = CLperf_prev_nonlumped;

        % info_mu_nonlumped contiene già la cronologia
        % del precedente controllore, caricata dal MAT.

        fprintf(['Nuovo controllore non-lumped scartato: ' ...
            'mantengo quello precedentemente salvato.\n']);

    end


    fprintf('Ordine K lumped     = %d\n', ...
        order(K_mu_scaled_lumped));

    fprintf('Ordine K non-lumped = %d\n', ...
        order(K_mu_scaled_nonlumped));

    fprintf('Performance selezionata lumped     = %.6f\n', ...
        CLperf_mu_lumped);

    fprintf('Performance selezionata non-lumped = %.6f\n', ...
        CLperf_mu_nonlumped);

else

    fprintf('\nSintesi diretta sul modello NON-LUMPED disattivata.\n');
    K_mu_scaled_nonlumped = [];
    CLperf_mu_nonlumped = NaN;
    info_mu_nonlumped = [];

end

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
%     J_alpha, l, omega_n, tau_d, Delta_Pade
%
% Pertanto anche la riduzione dell'ordine viene accettata
% solo se mantiene RP < 1 sul modello NON-LUMPED.

fprintf('\n============================================================\n');
fprintf('VERIFICA RP E RIDUZIONE SUL MODELLO NON-LUMPED\n');
fprintf('============================================================\n');

% Controllore full-order prodotto dalla sintesi lumped
K_mu_scaled_full = ss(K_mu_scaled_lumped);
N = order(K_mu_scaled_full);

fprintf('Ordine controllore full-order LUMPED = %d\n',N);
fprintf('RP della sintesi sul lumped           = %.6f\n',CLperf_mu_lumped);

%% RIDUZIONE PRELIMINARE CON MINREAL SU ENTRAMBI I CONTROLLORI
% Minreal elimina soltanto stati non controllabili/non osservabili o
% cancellazioni pole-zero; non cambia intenzionalmente la funzione di
% trasferimento ingresso-uscita. I candidati vengono validati separatamente
% sul modello incerto NON-LUMPED prima di qualsiasi eventuale selezione.

fprintf('\n============================================================\n');
fprintf('MINREAL SU ENTRAMBI I CONTROLLORI\n');
fprintf('============================================================\n');

K_mu_scaled_lumped_minreal = minreal(ss(K_mu_scaled_lumped));
N_lumped_full = order(K_mu_scaled_lumped);
N_lumped_minreal = order(K_mu_scaled_lumped_minreal);

fprintf('LUMPED: ordine originale = %d | dopo minreal = %d | stati eliminati = %d\n', ...
    N_lumped_full, N_lumped_minreal, N_lumped_full-N_lumped_minreal);

if RUN_MUSYN_NONLUMPED
    K_mu_scaled_nonlumped_minreal = minreal(ss(K_mu_scaled_nonlumped));
    N_nonlumped_full = order(K_mu_scaled_nonlumped);
    N_nonlumped_minreal = order(K_mu_scaled_nonlumped_minreal);

    fprintf('NON-LUMPED: ordine originale = %d | dopo minreal = %d | stati eliminati = %d\n', ...
        N_nonlumped_full, N_nonlumped_minreal, N_nonlumped_full-N_nonlumped_minreal);
else
    K_mu_scaled_nonlumped_minreal = [];
    N_nonlumped_full = NaN;
    N_nonlumped_minreal = NaN;
    fprintf('NON-LUMPED: minreal saltato perché la sintesi diretta è disattivata.\n');
end

% I due controllori minreal sono solo candidati di confronto: non sostituiscono
% automaticamente il controllore full-order selezionato più avanti.


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

%% CLOSED-LOOP NON-LUMPED CON K SINTETIZZATO NON-LUMPED (OPZIONALE)

CL_full_nonlumped_from_KN = [];

if RUN_MUSYN_NONLUMPED

    Knonlumped_full = ss(K_mu_scaled_nonlumped);

    Knonlumped_full.InputName  = {'e1','e2'};
    Knonlumped_full.OutputName = {'u1','u2'};

    CL_full_nonlumped_from_KN = connect( ...
        Gmu_nonlumped, ...
        Knonlumped_full, ...
        W_S, ...
        W_U, ...
        W_T, ...
        sum1, ...
        sum2, ...
        {'r1','r2'}, ...
        {'zS1','zS2','zU1','zU2','zT1','zT2'});

end

%% CLOSED-LOOP DEI CANDIDATI MINREAL SULLO STESSO MODELLO NON-LUMPED

K_lumped_minreal_test = ss(K_mu_scaled_lumped_minreal);
K_lumped_minreal_test.InputName  = {'e1','e2'};
K_lumped_minreal_test.OutputName = {'u1','u2'};

CL_minreal_lumped_nonlumped = connect( ...
    Gmu_nonlumped, ...
    K_lumped_minreal_test, ...
    W_S, W_U, W_T, sum1, sum2, ...
    {'r1','r2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2'});

CL_minreal_nonlumped_nonlumped = [];
if RUN_MUSYN_NONLUMPED
    K_nonlumped_minreal_test = ss(K_mu_scaled_nonlumped_minreal);
    K_nonlumped_minreal_test.InputName  = {'e1','e2'};
    K_nonlumped_minreal_test.OutputName = {'u1','u2'};

    CL_minreal_nonlumped_nonlumped = connect( ...
        Gmu_nonlumped, ...
        K_nonlumped_minreal_test, ...
        W_S, W_U, W_T, sum1, sum2, ...
        {'r1','r2'}, ...
        {'zS1','zS2','zU1','zU2','zT1','zT2'});
end


%% ------------------------------------------------------------------------
% 7.2 ROBUST PERFORMANCE FULL-ORDER NON-LUMPED
% -------------------------------------------------------------------------

fprintf('\nVerifica RP del controllore LUMPED sul modello NON-LUMPED...\n');

optsWC = wcOptions( ...
    'Display','off', ...
    'MussvOptions','a');

% UNICA verifica full-order indipendente mantenuta:
% il controllore sintetizzato sul lumped viene validato sul non-lumped.
[wc_full,wcu_full] = wcgain( ...
    CL_full_nonlumped, ...
    optsWC);

RP_full_nonlumped = wc_full.UpperBound;

fprintf('\nRP del K LUMPED sul NON-LUMPED:\n');
fprintf('  Lower bound = %.6f\n',wc_full.LowerBound);
fprintf('  Upper bound = %.6f\n',wc_full.UpperBound);
fprintf('  Frequenza critica = %.6f rad/s\n', ...
    wc_full.CriticalFrequency);

% La sintesi diretta non-lumped è già valutata da MUSYN.
% Non eseguiamo qui una seconda wcgain sul suo full-order.
fprintf('\nPerformance MUSYN del K direttamente NON-LUMPED:\n');

if RUN_MUSYN_NONLUMPED
    fprintf('  CLperf MUSYN = %.6f\n',CLperf_mu_nonlumped);
else
    fprintf('  Sintesi diretta disattivata.\n');
end

fprintf('\nWorst-case uncertainties del K LUMPED:\n');
disp(wcu_full);

%% WCGain DEL CONTROLLORE NON-LUMPED FULL-ORDER

if RUN_MUSYN_NONLUMPED

    fprintf('\n============================================================\n');
    fprintf('WCGain NON-LUMPED FULL-ORDER, PRIMA DI MINREAL\n');
    fprintf('============================================================\n');

    [wc_full_nonlumped_KN, wcu_full_nonlumped_KN] = ...
        wcgain(CL_full_nonlumped_from_KN, optsWC);

    RP_full_nonlumped_KN = wc_full_nonlumped_KN.UpperBound;

    fprintf('Ordine K non-lumped full-order = %d\n', ...
        order(K_mu_scaled_nonlumped));

    fprintf('RP lower = %.6f\n', ...
        wc_full_nonlumped_KN.LowerBound);

    fprintf('RP upper = %.6f\n', ...
        wc_full_nonlumped_KN.UpperBound);

    if RP_full_nonlumped_KN < 1
        fprintf('RP full-order NON-LUMPED: PASS\n');
    else
        fprintf('RP full-order NON-LUMPED: NON CERTIFICATA\n');
    end

end

%% VALIDAZIONE RP DEI CANDIDATI MINREAL
% Se minreal non ha eliminato stati, si riutilizza il risultato del full-order:
% la realizzazione ingresso-uscita è equivalente e si evita una wcgain inutile.

fprintf('\n============================================================\n');
fprintf('VALIDAZIONE RP DEI CANDIDATI MINREAL\n');
fprintf('============================================================\n');

if N_lumped_minreal < N_lumped_full
    fprintf('Verifica wcgain del controllore LUMPED dopo minreal...\n');
    [wc_minreal_lumped,wcu_minreal_lumped] = wcgain( ...
        CL_minreal_lumped_nonlumped, optsWC);
else
    fprintf('LUMPED: nessuna riduzione d''ordine con minreal; riuso wcgain full-order.\n');
    wc_minreal_lumped = wc_full;
    wcu_minreal_lumped = wcu_full;
end
RP_minreal_lumped = wc_minreal_lumped.UpperBound;

fprintf('LUMPED minreal: ordine %d -> %d | RP upper = %.6f\n', ...
    N_lumped_full, N_lumped_minreal, RP_minreal_lumped);

wc_minreal_nonlumped = [];
wcu_minreal_nonlumped = [];
RP_minreal_nonlumped = NaN;

if RUN_MUSYN_NONLUMPED

    fprintf('Verifica wcgain del controllore NON-LUMPED dopo minreal...\n');

    [wc_minreal_nonlumped,wcu_minreal_nonlumped] = wcgain( ...
        CL_minreal_nonlumped_nonlumped, ...
        optsWC);

    RP_minreal_nonlumped = wc_minreal_nonlumped.UpperBound;

    fprintf('NON-LUMPED minreal: ordine %d -> %d | RP lower = %.6f | RP upper = %.6f\n', ...
        N_nonlumped_full, ...
        N_nonlumped_minreal, ...
        wc_minreal_nonlumped.LowerBound, ...
        RP_minreal_nonlumped);
end

fprintf('\n============================================================\n');
fprintf('VALIDAZIONE SUL MODELLO NON-LUMPED\n');
fprintf('============================================================\n');

fprintf('K sintetizzato LUMPED:\n');
fprintf('  RP lower = %.6f\n',wc_full.LowerBound);
fprintf('  RP upper = %.6f\n',wc_full.UpperBound);

if RUN_MUSYN_NONLUMPED
    fprintf('\nK sintetizzato direttamente NON-LUMPED:\n');
    fprintf('  CLperf MUSYN = %.6f\n',CLperf_mu_nonlumped);

    fprintf('\nK NON-LUMPED dopo minreal:\n');
    fprintf('  RP lower = %.6f\n',wc_minreal_nonlumped.LowerBound);
    fprintf('  RP upper = %.6f\n',RP_minreal_nonlumped);

    fprintf('\nWorst-case uncertainties del candidato non-lumped dopo minreal:\n');
    disp(wcu_minreal_nonlumped);
else
    fprintf('\nVerifica del controllore sintetizzato direttamente sul non-lumped: saltata.\n');
end


fprintf('\n============================================================\n');
fprintf('WORST-CASE UNCERTAINTIES\n');
fprintf('============================================================\n');

fprintf('Worst-case J_alpha = %.12g\n', wcu_full.J_alpha);
fprintf('Worst-case l       = %.12g\n', wcu_full.l);
fprintf('\nWorst-case uncertainties:\n');
disp(wcu_full);

RP_full_nonlumped = wc_full.UpperBound;

fprintf('RP full-order NON-LUMPED lower = %.6f\n', ...
    wc_full.LowerBound);

fprintf('RP full-order NON-LUMPED upper = %.6f\n', ...
    wc_full.UpperBound);

fprintf('Frequenza critica = %.6f rad/s\n', ...
    wc_full.CriticalFrequency);



fprintf('RP full-order NON-LUMPED = %.6f\n', ...
    RP_full_nonlumped);

%% ========================================================================
% SALVATAGGIO DEI MIGLIORI CONTROLLORI PER LE PROSSIME RIPRESE
% ========================================================================

K_mu_scaled_lumped_resume = ss(K_mu_scaled_lumped);
CLperf_mu_lumped_resume = CLperf_mu_lumped;

if RUN_MUSYN_NONLUMPED

    K_mu_scaled_nonlumped_resume = ss(K_mu_scaled_nonlumped);
    CLperf_mu_nonlumped_resume = CLperf_mu_nonlumped;

elseif isfield(resumeData,'K_mu_scaled_nonlumped_resume')

    % Se la sintesi non-lumped è disattivata, conserva il precedente risultato.
    K_mu_scaled_nonlumped_resume = ...
        resumeData.K_mu_scaled_nonlumped_resume;

    CLperf_mu_nonlumped_resume = ...
        resumeData.CLperf_mu_nonlumped_resume;

else

    K_mu_scaled_nonlumped_resume = ...
        resumeData.K_mu_scaled_nonlumped;

    CLperf_mu_nonlumped_resume = ...
        resumeData.CLperf_mu_nonlumped;

end

save('MU_comparison_lumped_nonlumped.mat', ...
    'K_mu_scaled_lumped_resume', ...
    'CLperf_mu_lumped_resume', ...
    'K_mu_scaled_nonlumped_resume', ...
    'CLperf_mu_nonlumped_resume', ...
    'K_mu_scaled_lumped', ...
    'CLperf_mu_lumped', ...
    'info_mu_lumped', ...
    'K_mu_scaled_nonlumped', ...
    'CLperf_mu_nonlumped', ...
    'info_mu_nonlumped', ...
    'wc_full', ...
    'wcu_full', ...
    'RP_full_nonlumped', ...
    'K_mu_scaled_lumped_minreal', ...
    'K_mu_scaled_nonlumped_minreal', ...
    'N_lumped_full', 'N_lumped_minreal', ...
    'N_nonlumped_full', 'N_nonlumped_minreal', ...
    'wc_minreal_lumped', 'wcu_minreal_lumped', 'RP_minreal_lumped', ...
    'wc_minreal_nonlumped', 'wcu_minreal_nonlumped', 'RP_minreal_nonlumped');

fprintf('Confronto lumped/non-lumped salvato.\n');

if RP_full_nonlumped < 1
    fprintf('K lumped verificato sul non-lumped: RP < 1 [PASS]\n');
else
    warning(['K sintetizzato sul lumped NON soddisfa RP < 1 ' ...
             'quando viene verificato sul modello non-lumped.']);
end

if RUN_MUSYN_NONLUMPED
    if RP_minreal_nonlumped < 1
        fprintf('K non-lumped dopo minreal: RP < 1 [PASS]\n');
    else
        warning(['Il candidato non-lumped dopo minreal ' ...
            'NON soddisfa RP < 1.']);
    end
end

fprintf('\n============================================================\n');
fprintf('RIDUZIONE OTTIMIZZATA DELL''ORDINE\n');
fprintf('============================================================\n');

K_mu_scaled_reduced = [];
RP_reduced = NaN;
order_reduced = NaN;
Kred = [];

%% ------------------------------------------------------------------------
% 7.3 GENERAZIONE DEI CANDIDATI PER LA RIDUZIONE (OPZIONALE)
% -------------------------------------------------------------------------

if RUN_ORDER_REDUCTION

    orders_to_test = 8:(N-1);

    if ~isempty(orders_to_test)

        % Calcola i controllori ridotti per tutti gli ordini richiesti
        Kred = balred(K_mu_scaled_full, orders_to_test);

        % Verifica ciascun controllore ridotto
        for idx = 1:numel(orders_to_test)

            k = orders_to_test(idx);
            Ktest = ss(Kred(:,:,idx));

            Ktest.InputName  = {'e1','e2'};
            Ktest.OutputName = {'u1','u2'};

            % Closed-loop incerto NON-LUMPED
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

            % 1) Controllo NOMINALE
            CLnom_test = getNominal(CLtest);

            if ~isstable(CLnom_test)
                fprintf('Ordine %3d --> NP = Inf --> skip\n',k);
                continue;
            end

            NP_test = hinfnorm(minreal(CLnom_test,1e-5));

            % Se NP >= 1, RP non può essere < 1
            if NP_test >= 1
                fprintf('Ordine %3d --> NP = %.6f --> skip WCgain\n', ...
                    k,NP_test);
                continue;
            end

            % 2) Solo se NP < 1 facciamo wcgain
            [wc_test,~] = wcgain(CLtest,optsWC);
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
    else
        fprintf('Nessun ordine ridotto da testare (N = %d).\n',N);
    end

else
    fprintf('\nRicerca del controllore ridotto disattivata.\n');
end

%% ------------------------------------------------------------------------
% 7.4 SELEZIONE CONTROLLORE FINALE
% -------------------------------------------------------------------------

if ~RUN_ORDER_REDUCTION

    % Nessuna riduzione richiesta: mantieni il controllore full-order.
    K_mu_scaled = K_mu_scaled_full;
    RP_final = RP_full_nonlumped;
    order_final = N;

    fprintf('Riduzione disattivata: mantenuto il controllore full-order.\n');
    fprintf('RP full-order sul NON-LUMPED = %.6f\n',RP_final);

elseif isempty(K_mu_scaled_reduced)

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

    fprintf('Ordine controllore finale = %d\n',order_final);
    fprintf('RP finale NON-LUMPED       = %.6f\n',RP_final);
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
% 10. CONFRONTO DELLE PRESTAZIONI NOMINALI SUL NON-LUMPED
% ========================================================================

G_nom_nonlumped = G_unc_mu_scaled.NominalValue;
I2 = eye(2);

%% 10.1 Controllore sintetizzato sul LUMPED, ordine pieno
% Si usa la stessa interconnessione connect già validata con wcgain.
% CL_full_nonlumped contiene il K full-order sintetizzato sul lumped,
% applicato al modello parametrico non-lumped.

CLnom_KL = getNominal(CL_full_nonlumped);
CLnom_KL = minreal(CLnom_KL,1e-5);

if isstable(CLnom_KL)
    gammaNom_lumped_synth = hinfnorm(CLnom_KL);
else
    gammaNom_lumped_synth = Inf;
    warning('Il closed-loop nominale del K lumped sul non-lumped è instabile.');
end

%% 10.2 Controllore sintetizzato direttamente sul NON-LUMPED (OPZIONALE)

gammaNom_nonlumped_synth = NaN;

if RUN_MUSYN_NONLUMPED

    CLnom_KN = getNominal(CL_full_nonlumped_from_KN);
    CLnom_KN = minreal(CLnom_KN,1e-5);

    if isstable(CLnom_KN)
        gammaNom_nonlumped_synth = hinfnorm(CLnom_KN);
    else
        gammaNom_nonlumped_synth = Inf;
        warning('Il closed-loop nominale del K non-lumped è instabile.');
    end

end

%% 10.3 Controllore FINALE selezionato, eventualmente ridotto
% Ricostruisce il closed-loop pesato con il controllore effettivamente
% selezionato in sezione 7.4. Si conserva anche Tnom_nonlumped per la
% verifica separata dei poli nella sezione 11.

Lnom_nonlumped = G_nom_nonlumped * K_mu_scaled;
Snom_nonlumped = feedback(I2,Lnom_nonlumped);
Tnom_nonlumped = feedback(Lnom_nonlumped,I2);

Kfinal_test = ss(K_mu_scaled);
Kfinal_test.InputName  = {'e1','e2'};
Kfinal_test.OutputName = {'u1','u2'};

CL_final_nonlumped = connect( ...
    Gmu_nonlumped, ...
    Kfinal_test, ...
    W_S, ...
    W_U, ...
    W_T, ...
    sum1, ...
    sum2, ...
    {'r1','r2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2'});

CLnom_nonlumped = getNominal(CL_final_nonlumped);
CLnom_nonlumped = minreal(CLnom_nonlumped,1e-5);

if isstable(CLnom_nonlumped)
    gammaNom_nonlumped = hinfnorm(CLnom_nonlumped);
else
    gammaNom_nonlumped = Inf;
    warning('Il closed-loop nominale del controllore finale è instabile.');
end

%% Stampa confronto

fprintf('\n============================================================\n');
fprintf('CONFRONTO PRESTAZIONI NOMINALI SULLO STESSO NON-LUMPED\n');
fprintf('============================================================\n');

fprintf('Gamma K sintetizzato LUMPED, ordine pieno = %.6f\n', ...
    gammaNom_lumped_synth);

if RUN_MUSYN_NONLUMPED
    fprintf('Gamma K sintetizzato NON-LUMPED          = %.6f\n', ...
        gammaNom_nonlumped_synth);
else
    fprintf('Gamma K sintetizzato NON-LUMPED          = non calcolato\n');
end

fprintf('Gamma K FINALE selezionato               = %.6f\n', ...
    gammaNom_nonlumped);

% La NP finale viene valutata sul modello NON-LUMPED nella sezione 10.3.

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
% CONTROLLORI FISICI DEI DUE RAMI DI SINTESI
% ========================================================================

% Controllore sintetizzato sul modello LUMPED, ordine pieno
K_mu_lumped_phys = Du * K_mu_scaled_lumped * Dy_inv;
K_mu_lumped_phys.InputName  = {'e1','e2'};
K_mu_lumped_phys.OutputName = {'u1','u2'};

% Candidato LUMPED ottenuto con minreal (non selezionato automaticamente)
K_mu_lumped_minreal_phys = Du * K_mu_scaled_lumped_minreal * Dy_inv;
K_mu_lumped_minreal_phys.InputName  = {'e1','e2'};
K_mu_lumped_minreal_phys.OutputName = {'u1','u2'};

% Controllore sintetizzato direttamente sul NON-LUMPED (opzionale)
if RUN_MUSYN_NONLUMPED

    K_mu_nonlumped_phys = Du * K_mu_scaled_nonlumped * Dy_inv;
    K_mu_nonlumped_phys.InputName  = {'e1','e2'};
    K_mu_nonlumped_phys.OutputName = {'u1','u2'};

    K_mu_nonlumped_minreal_phys = Du * K_mu_scaled_nonlumped_minreal * Dy_inv;
    K_mu_nonlumped_minreal_phys.InputName  = {'e1','e2'};
    K_mu_nonlumped_minreal_phys.OutputName = {'u1','u2'};

    fprintf('\nControllori conservati entrambi:\n');
    fprintf('K_mu_scaled_lumped:     ordine %d\n', ...
        order(K_mu_scaled_lumped));
    fprintf('K_mu_scaled_nonlumped:  ordine %d\n', ...
        order(K_mu_scaled_nonlumped));

else

    K_mu_nonlumped_phys = [];
    K_mu_nonlumped_minreal_phys = [];
    fprintf('\nControllore fisico non-lumped non creato: sintesi disattivata.\n');
    fprintf('K_mu_scaled_lumped:     ordine %d\n', ...
        order(K_mu_scaled_lumped));

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
    'P_mu',...
    'K_mu_scaled_lumped', ...
    'CLperf_mu_lumped', ...
    'info_mu_lumped', ...
    'K_mu_scaled_nonlumped', ...
    'CLperf_mu_nonlumped', ...
    'info_mu_nonlumped', ...
    'K_mu_scaled_lumped_resume', ...
    'CLperf_mu_lumped_resume', ...
    'K_mu_scaled_nonlumped_resume', ...
    'CLperf_mu_nonlumped_resume', ...
    'P_mu_nonlumped', ...
    'K_mu_lumped_phys', ...
    'K_mu_nonlumped_phys', ...
    'K_mu_scaled_lumped_minreal', ...
    'K_mu_scaled_nonlumped_minreal', ...
    'K_mu_lumped_minreal_phys', ...
    'K_mu_nonlumped_minreal_phys', ...
    'N_lumped_full', 'N_lumped_minreal', ...
    'N_nonlumped_full', 'N_nonlumped_minreal', ...
    'RP_minreal_lumped', 'RP_minreal_nonlumped', ...
    'wc_minreal_lumped', 'wcu_minreal_lumped', ...
    'wc_minreal_nonlumped', 'wcu_minreal_nonlumped', ...
    'gammaNom_lumped_synth', ...
    'gammaNom_nonlumped_synth', ...
    'wc_full', ...
    'wcu_full');

fprintf('\n============================================================\n');
fprintf('SALVATAGGIO COMPLETATO\n');
fprintf('============================================================\n');

fprintf('File "MU_controller.mat" creato correttamente.\n');


%% ========================================================================
% 6.2 VERIFICHE NS, NP, RS, RP DEL CONTROLL0RE WARM-START MIXSYN
% ========================================================================
%
% Il controllore K_direct_mixsyn è sintetizzato sul modello
% non-lumped semplificato, ma viene validato sul modello originale:
%
%       G_unc_mu_scaled
%
% Le verifiche sono effettuate nelle coordinate scalate,
% utilizzando WS, WU e WT.
%
% RP viene valutata sulla stessa griglia omegaWeights
% utilizzata nel confronto dei controllori.

fprintf('\n============================================================\n');
fprintf('VALIDAZIONE K WARM-START MIXSYN SUL NON-LUMPED ORIGINALE\n');
fprintf('============================================================\n');

%% 1. PLANT E CONTROLL0RE

Gval = G_unc_mu_scaled;
Gval.InputName = {'u1','u2'};
Gval.OutputName = {'y1','y2'};

Kval = ss(K_direct_mixsyn);
Kval.InputName = {'e1','e2'};
Kval.OutputName = {'u1','u2'};

I2 = eye(2);

%% 2. PESI

W_S_val = ss(WS);
W_S_val.InputName = {'e1','e2'};
W_S_val.OutputName = {'zS1','zS2'};

W_U_val = ss(WU);
W_U_val.InputName = {'u1','u2'};
W_U_val.OutputName = {'zU1','zU2'};

W_T_val = ss(WT);
W_T_val.InputName = {'y1','y2'};
W_T_val.OutputName = {'zT1','zT2'};

%% 3. INTERCONNESSIONE PESATA INCERTA

Sum1_val = sumblk('e1 = r1 - y1');
Sum2_val = sumblk('e2 = r2 - y2');

CL_unc_direct_mixsyn = connect( ...
    Gval, Kval, ...
    W_S_val, W_U_val, W_T_val, ...
    Sum1_val, Sum2_val, ...
    {'r1','r2'}, ...
    {'zS1','zS2','zU1','zU2','zT1','zT2'});

CL_nom_direct_mixsyn = getNominal(CL_unc_direct_mixsyn);

%% 4. NS - NOMINAL STABILITY

Tnom_direct_mixsyn = ...
    feedback(Gval.NominalValue*K_direct_mixsyn,I2);

NS_direct_mixsyn = isstable(Tnom_direct_mixsyn);

fprintf('\n--- NOMINAL STABILITY (NS) ---\n');
fprintf('NS = %d\n',NS_direct_mixsyn);

if NS_direct_mixsyn
    fprintf('NS PASS: sistema nominalmente stabile.\n');
else
    fprintf('NS FAIL: sistema nominalmente instabile.\n');
end

%% 5. NP - NOMINAL PERFORMANCE
% Calcolo senza minreal, sul sistema pesato nominale.

if NS_direct_mixsyn

    [NP_direct_mixsyn,fNP_direct_mixsyn] = ...
        hinfnorm(CL_nom_direct_mixsyn);

else

    NP_direct_mixsyn = Inf;
    fNP_direct_mixsyn = NaN;

end

fprintf('\n--- NOMINAL PERFORMANCE (NP) ---\n');
fprintf('Gamma NP = %.9f\n',NP_direct_mixsyn);
fprintf('Frequenza di picco = %.6f rad/s\n', ...
    fNP_direct_mixsyn);

if NP_direct_mixsyn < 1
    fprintf('NP PASS: gamma nominale < 1.\n');
else
    fprintf('NP FAIL: gamma nominale >= 1.\n');
end

%% 6. RS - ROBUST STABILITY

fprintf('\n--- ROBUST STABILITY (RS) ---\n');

optsRS_direct_mixsyn = ...
    robuststabOptions('Display','off');

Tunc_direct_mixsyn = ...
    feedback(Gval*K_direct_mixsyn,I2);

[stabmarg_direct_mixsyn, ...
    destabunc_direct_mixsyn, ...
    reportRS_direct_mixsyn, ...
    infoRS_direct_mixsyn] = ...
    robuststab(Tunc_direct_mixsyn,optsRS_direct_mixsyn);

RS_lower_direct_mixsyn = ...
    stabmarg_direct_mixsyn.LowerBound;

RS_upper_direct_mixsyn = ...
    stabmarg_direct_mixsyn.UpperBound;

fprintf('RS lower bound = %.9f\n',RS_lower_direct_mixsyn);
fprintf('RS upper bound = %.9f\n',RS_upper_direct_mixsyn);

if RS_lower_direct_mixsyn > 1
    fprintf('RS PASS: stabilita robusta certificata.\n');
elseif RS_upper_direct_mixsyn < 1
    fprintf('RS FAIL: stabilita robusta violata.\n');
else
    fprintf('RS INCONCLUSIVA: bound comprendono 1.\n');
end

%% 7. RP - ROBUST PERFORMANCE
% Stessa griglia omegaWeights usata nel confronto precedente.

fprintf('\n--- ROBUST PERFORMANCE (RP) ---\n');

optsRP_direct_mixsyn = wcOptions( ...
    'Display','off', ...
    'MussvOptions','a');

omegaRP_direct_mixsyn = omegaWeights;

[wc_direct_mixsyn,wcu_direct_mixsyn] = ...
    wcgain( ...
        ufrd(CL_unc_direct_mixsyn,omegaRP_direct_mixsyn), ...
        optsRP_direct_mixsyn);

RP_lower_direct_mixsyn = ...
    wc_direct_mixsyn.LowerBound;

RP_upper_direct_mixsyn = ...
    wc_direct_mixsyn.UpperBound;

fprintf('RP lower bound = %.9f\n',RP_lower_direct_mixsyn);
fprintf('RP upper bound = %.9f\n',RP_upper_direct_mixsyn);

if RP_upper_direct_mixsyn < 1
    fprintf('RP PASS: prestazione robusta certificata.\n');
elseif RP_lower_direct_mixsyn > 1
    fprintf('RP FAIL: prestazione robusta violata.\n');
else
    fprintf('RP INCONCLUSIVA: bound comprendono 1.\n');
end

%% 8. RIEPILOGO

fprintf('\n============================================================\n');
fprintf('RIEPILOGO NS, NP, RS, RP - WARM-START MIXSYN\n');
fprintf('============================================================\n');

fprintf('Ordine K = %d\n',order(K_direct_mixsyn));
fprintf('CLperf MUSYN = %.9f\n',CLperf_direct_mixsyn);

fprintf('NS = %d\n',NS_direct_mixsyn);
fprintf('NP = %.9f | PASS = %d\n', ...
    NP_direct_mixsyn,NP_direct_mixsyn < 1);

fprintf('RS = [%.9f, %.9f] | PASS = %d\n', ...
    RS_lower_direct_mixsyn, ...
    RS_upper_direct_mixsyn, ...
    RS_lower_direct_mixsyn > 1);

fprintf('RP = [%.9f, %.9f] | PASS = %d\n', ...
    RP_lower_direct_mixsyn, ...
    RP_upper_direct_mixsyn, ...
    RP_upper_direct_mixsyn < 1);

%% 9. SALVATAGGIO DEI RISULTATI SENZA SOVRASCRIVERE
% -append conserva K, cronologia D-K e opzioni già salvati.

save('MU_DIRECT_FROM_MIXSYN_INIT.mat', ...
    'NS_direct_mixsyn', ...
    'NP_direct_mixsyn', ...
    'fNP_direct_mixsyn', ...
    'stabmarg_direct_mixsyn', ...
    'destabunc_direct_mixsyn', ...
    'reportRS_direct_mixsyn', ...
    'infoRS_direct_mixsyn', ...
    'RS_lower_direct_mixsyn', ...
    'RS_upper_direct_mixsyn', ...
    'wc_direct_mixsyn', ...
    'wcu_direct_mixsyn', ...
    'RP_lower_direct_mixsyn', ...
    'RP_upper_direct_mixsyn', ...
    'CL_unc_direct_mixsyn', ...
    'CL_nom_direct_mixsyn', ...
    '-append');

fprintf('\nVerifiche salvate nel MAT senza eliminare i dati precedenti.\n');
