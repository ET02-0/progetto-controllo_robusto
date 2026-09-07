%% ========================================================================
% MU_ANALYSIS_COMPARISON
%
% Confronto H-infinity vs mu-synthesis sul modello incerto del progetto.
% Adattato ai modelli e ai controllori del progetto corrente.
% ========================================================================

close all;
clc;

%% ========================================================================
% 0. CARICAMENTO DATI
% ========================================================================

load('HINF_workspace.mat');
load('HINF_controllers.mat');
load('ACTUATOR_LUMPED.mat');
load('MU_controller.mat');

I2 = eye(2);

%% ========================================================================
% 1. CONTROLLORI DA CONFRONTARE
% ========================================================================

controllers = {
    K_mix_scaled
    K_hinfsyn_scaled
    K_pidcomp_scaled
    K_mu_scaled
};

controllerNames = {
    'mixsyn'
    'hinfsyn'
    'PID+comp'
    'mu-synthesis'
};

Nc = numel(controllers);

%% ========================================================================
% 2. PIANTE PER L'ANALISI
% ========================================================================

% Plant nominale per NS e NP
Gnom = G_scaled;

% Plant incerto lumped per RS, RP e analisi mu
Gunc = G_uncertain_lumped_scaled;

Gunc.InputName = {
    'u1'
    'u2'
};

Gunc.OutputName = {
    'y1'
    'y2'
};
%% ========================================================================
% 3. GRIGLIA FREQUENZIALE PER ANALISI MU
% ========================================================================

omegaAnalysis = omegaWeights

fprintf('\n============================================================\n');
fprintf('MU ANALYSIS COMPARISON\n');
fprintf('============================================================\n');

fprintf('Numero frequenze = %d\n', ...
    numel(omegaAnalysis));

fprintf('Intervallo = %.3e - %.3e rad/s\n', ...
    min(omegaAnalysis), ...
    max(omegaAnalysis));

%% ========================================================================
% 4. OPZIONI WCGain
% ========================================================================

wcOptsFast = wcOptions;

wcOptsAccurate = wcOptions( ...
    'MussvOptions','a');

%% ========================================================================
% 5. PREALLOCAZIONE
% ========================================================================

NominalGamma = nan(Nc,1);

RS_Lower = nan(Nc,1);
RS_Upper = nan(Nc,1);

RP_Lower = nan(Nc,1);
RP_Upper = nan(Nc,1);

MuRS_Lower = nan(Nc,1);
MuRS_Upper = nan(Nc,1);

MuRP_Lower = nan(Nc,1);
MuRP_Upper = nan(Nc,1);

NS = false(Nc,1);
NP = false(Nc,1);
RS = false(Nc,1);
RP = false(Nc,1);

%% ========================================================================
% 6. ANALISI DEI CONTROLLORI
% ========================================================================

for k = 1:Nc

    K = controllers{k};

    fprintf('\n============================================================\n');
    fprintf('Controller %d/%d: %s\n', ...
        k,Nc,controllerNames{k});
    fprintf('============================================================\n');

    %% ====================================================================
    % 6.1 ANALISI NOMINALE
    % =====================================================================

    fprintf('\nNominal analysis...\n');

    Lnom = Gnom * K;

    Snom = feedback(I2,Lnom);
    Tnom = feedback(Lnom,I2);

    Wnom = [
        WS*Snom
        WU*K*Snom
        WT*Tnom
    ];

    NominalGamma(k) = ...
        hinfnorm( ...
            minreal(Wnom,1e-7));

    NS(k) = ...
        all(real(pole(Tnom)) < 0);

    NP(k) = ...
        NominalGamma(k) < 1;

    fprintf('Nominal Gamma = %.6f\n', ...
        NominalGamma(k));

    fprintf('NS = %d\n',NS(k));
    fprintf('NP = %d\n',NP(k));

    %% ====================================================================
    % 6.2 SISTEMA INCERTO
    % =====================================================================

    Lunc = Gunc * K;

    Sunc = feedback(I2,Lunc);
    Tunc = feedback(Lunc,I2);

    Wunc = [
        WS*Sunc
        WU*K*Sunc
        WT*Tunc
    ];

    %% ====================================================================
    % 6.3 ROBUST STABILITY - ROBUSTSTAB
    % =====================================================================

    fprintf('\nRunning robstab...\n');

    tic;

    stabMargin = ...
        robstab(Tunc);

    elapsedRobstab = toc;

    RS_Lower(k) = ...
        stabMargin.LowerBound;

    RS_Upper(k) = ...
        stabMargin.UpperBound;

    RS(k) = ...
        RS_Lower(k) > 1;

    fprintf('robstab completed in %.2f s.\n', ...
        elapsedRobstab);

    fprintf('RS Lower = %.6f\n', ...
        RS_Lower(k));

    fprintf('RS Upper = %.6f\n', ...
        RS_Upper(k));

    %% ====================================================================
    % 6.4 ROBUST PERFORMANCE - WCGain
    % =====================================================================

    fprintf('\nRunning wcgain...\n');

    tic;

    try

        wcGain = ...
            wcgain( ...
                Wunc, ...
                wcOptsFast);

        fprintf('Standard wcgain successful.\n');

    catch ME

        if contains( ...
                ME.message, ...
                'Invalid MU upper bound', ...
                'IgnoreCase',true)

            fprintf('\n');
            fprintf('Standard wcgain failed.\n');
            fprintf('Retrying with MussvOptions = ''a''...\n');

            wcGain = ...
                wcgain( ...
                    Wunc, ...
                    wcOptsAccurate);

            fprintf('Accurate wcgain successful.\n');

        else

            rethrow(ME);

        end

    end

    elapsedWcgain = toc;

    fprintf('wcgain completed in %.2f s.\n', ...
        elapsedWcgain);

    RP_Lower(k) = ...
        wcGain.LowerBound;

    RP_Upper(k) = ...
        wcGain.UpperBound;

    RP(k) = ...
        RP_Upper(k) < 1;

    fprintf('RP Lower = %.6f\n', ...
        RP_Lower(k));

    fprintf('RP Upper = %.6f\n', ...
        RP_Upper(k));

    %% ====================================================================
    % 6.5 DECOMPOSIZIONE LFT
    % =====================================================================

    fprintf('\nExtracting LFT data...\n');

    [Mdelta,Delta,BlockStructure] = ...
        lftdata(Wunc);

    %% ====================================================================
    % 6.6 MU ESPLICITA - ROBUST STABILITY
    % =====================================================================

    szDelta = size(Delta);

    M11 = ...
        Mdelta( ...
            1:szDelta(2), ...
            1:szDelta(1));

    fprintf( ...
        'Running explicit mu RS on %d frequencies...\n', ...
        numel(omegaAnalysis));

    tic;

    muRS = ...
        mussv( ...
            frd(M11,omegaAnalysis), ...
            BlockStructure, ...
            's');

    elapsedMuRS = toc;

    fprintf('mu RS completed in %.2f s.\n', ...
        elapsedMuRS);

    muRSup = squeeze( ...
        muRS.ResponseData(1,1,:));

    muRSlo = squeeze( ...
        muRS.ResponseData(1,2,:));

    MuRS_Upper(k) = ...
        max(muRSup);

    MuRS_Lower(k) = ...
        max(muRSlo);

    %% ====================================================================
    % 6.7 MU ESPLICITA - ROBUST PERFORMANCE
    % =====================================================================

    nExogenous = ...
        size(Wunc,2);

    nPerformance = ...
        size(Wunc,1);

    BlockStructureRP = ...
        BlockStructure;

    perfBlock = ...
        BlockStructure(1);

    perfBlock.Name = ...
        'DeltaPerf';

    perfBlock.Size = ...
        [nExogenous nPerformance];

    perfBlock.Type = ...
        'ultidyn';

    perfBlock.Occurrences = ...
        1;

    perfBlock.Simplify = ...
        0;

    BlockStructureRP(end+1,1) = ...
        perfBlock;

    fprintf( ...
        'Running explicit mu RP on %d frequencies...\n', ...
        numel(omegaAnalysis));

    tic;

    muRP = ...
        mussv( ...
            frd(Mdelta,omegaAnalysis), ...
            BlockStructureRP, ...
            's');

    elapsedMuRP = toc;

    fprintf('mu RP completed in %.2f s.\n', ...
        elapsedMuRP);

    muRPup = squeeze( ...
        muRP.ResponseData(1,1,:));

    muRPlo = squeeze( ...
        muRP.ResponseData(1,2,:));

    MuRP_Upper(k) = ...
        max(muRPup);

    MuRP_Lower(k) = ...
        max(muRPlo);

    %% ====================================================================
    % 6.8 INFORMAZIONI AGGIUNTIVE NOMINALI
    % =====================================================================

    [gWS_nom,wWS_nom] = ...
        hinfnorm( ...
            minreal(WS*Snom,1e-7));

    [gWU_nom,wWU_nom] = ...
        hinfnorm( ...
            minreal(WU*K*Snom,1e-7));

    [gWT_nom,wWT_nom] = ...
        hinfnorm( ...
            minreal(WT*Tnom,1e-7));

    fprintf('\nNominal weighted channels:\n');

    fprintf('WS*S  = %.6f at %.6f rad/s\n', ...
        gWS_nom,wWS_nom);

    fprintf('WU*KS = %.6f at %.6f rad/s\n', ...
        gWU_nom,wWU_nom);

    fprintf('WT*T  = %.6f at %.6f rad/s\n', ...
        gWT_nom,wWT_nom);

    %% ====================================================================
    % 6.9 RISULTATO CONTROLLER
    % =====================================================================

    fprintf('\nController %s completed.\n', ...
        controllerNames{k});

    fprintf('Nominal gamma : %.6f\n', ...
        NominalGamma(k));

    fprintf('RS margin     : [%.6f, %.6f]\n', ...
        RS_Lower(k), ...
        RS_Upper(k));

    fprintf('WC gain       : [%.6f, %.6f]\n', ...
        RP_Lower(k), ...
        RP_Upper(k));

    fprintf('mu RS         : [%.6f, %.6f]\n', ...
        MuRS_Lower(k), ...
        MuRS_Upper(k));

    fprintf('mu RP         : [%.6f, %.6f]\n', ...
        MuRP_Lower(k), ...
        MuRP_Upper(k));

end

%% ========================================================================
% 7. TABELLA RISULTATI
% ========================================================================

resultsMU = table( ...
    controllerNames, ...
    NominalGamma, ...
    RS_Lower, ...
    RS_Upper, ...
    RP_Lower, ...
    RP_Upper, ...
    MuRS_Lower, ...
    MuRS_Upper, ...
    MuRP_Lower, ...
    MuRP_Upper, ...
    NS, ...
    NP, ...
    RS, ...
    RP, ...
    'VariableNames',{ ...
        'Controller'
        'NominalGamma'
        'RS_Lower'
        'RS_Upper'
        'RP_Lower'
        'RP_Upper'
        'MuRS_Lower'
        'MuRS_Upper'
        'MuRP_Lower'
        'MuRP_Upper'
        'NS'
        'NP'
        'RS'
        'RP'
    });

%% ========================================================================
% 8. STAMPA RISULTATI
% ========================================================================

fprintf('\n============================================================\n');
fprintf('CONFRONTO HINF vs MU-SYNTHESIS\n');
fprintf('============================================================\n');

disp(resultsMU);

%% ========================================================================
% 9. SALVATAGGIO
% ========================================================================

writetable( ...
    resultsMU, ...
    'MU_vs_HINF_results.csv');

fprintf('\nRisultati salvati in MU_vs_HINF_results.csv\n');