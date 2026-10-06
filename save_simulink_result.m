%% ========================================================================
% SAVE_SIMULINK_RESULT
% ========================================================================
% Salva i risultati To Workspace dell'ULTIMA simulazione Simulink.
%
% WORKFLOW:
%   1) Eseguire MANUALMENTE Simulink con il controllore desiderato.
%   2) Impostare SOLO queste due righe:
%        simController = 'mu';
%   3) Eseguire QUESTO script.
%   4) Ripetere per gli altri controllori e/o per linear/nonlinear.
%   5) Alla fine eseguire UNA SOLA VOLTA MU_ANALYSIS_COMPARISON.m.
%
% Questo file NON esegue Simulink e NON crea i grafici comparativi.
% Serve solo a congelare ogni run manuale in uno snapshot MAT separato.
%
% To Workspace attesi:
%   MU:
%       alpha_mu beta_mu alpha_meas_mu beta_meas_mu u_cmd_mu u_sat_mu
%   HINF (mixsyn, hinfsyn, PID+comp):
%       alpha_hinf beta_hinf alpha_meas_hinf beta_meas_hinf
%       u_cmd_hinf u_sat_hinf
% ========================================================================


%% 1. CONFIGURAZIONE RUN CORRENTE
simController = 'mixsyn';          % <-- mu | mixsyn | hinfsyn | pidcomp


simController = lower(strtrim(simController));

validControllers = {'mu','mixsyn','hinfsyn','pidcomp'};

if ~ismember(simController,validControllers)
    error('simController non valido. Usa: mu | mixsyn | hinfsyn | pidcomp');
end

%% 2. CARTELLA RISULTATI
resultDir = fullfile(pwd,'simulink_results');
if ~exist(resultDir,'dir')
    mkdir(resultDir);
end

%% 3. NOMI TO WORKSPACE

if strcmp(simController,'mu')

    % LINEARE
    names.alpha     = 'alpha_mu';
    names.beta      = 'beta_mu';
    names.alphaMeas = 'alpha_meas_mu';
    names.betaMeas  = 'beta_meas_mu';
    names.uCmd      = 'u_cmd_mu';
    names.uSat      = 'u_sat_mu';

    % NON LINEARE
    names.alphaNL     = 'alpha_mu_nl';
    names.betaNL      = 'beta_mu_nl';
    names.alphaMeasNL = 'alpha_meas_mu_nl';
    names.betaMeasNL  = 'beta_meas_mu_nl';
    names.uCmdNL      = 'u_cmd_mu_nl';
    names.uSatNL      = 'u_sat_mu_nl';

else

    % LINEARE
    names.alpha     = 'alpha_hinf';
    names.beta      = 'beta_hinf';
    names.alphaMeas = 'alpha_meas_hinf';
    names.betaMeas  = 'beta_meas_hinf';
    names.uCmd      = 'u_cmd_hinf';
    names.uSat      = 'u_sat_hinf';

    % NON LINEARE
    names.alphaNL     = 'alpha_hinf_nl';
    names.betaNL      = 'beta_hinf_nl';
    names.alphaMeasNL = 'alpha_meas_hinf_nl';
    names.betaMeasNL  = 'beta_meas_hinf_nl';
    names.uCmdNL      = 'u_cmd_hinf_nl';
    names.uSatNL      = 'u_sat_hinf_nl';

end
%% 4. LETTURA SEGNALI DA out

if ~exist('out','var')
    error(['La variabile "out" non esiste nel Workspace. ', ...
           'Esegui prima la simulazione Simulink.']);
end

%% ------------------------------------------------------------
% LINEARE
% ------------------------------------------------------------

requiredLinear = { ...
    names.alpha, ...
    names.beta, ...
    names.uCmd, ...
    names.uSat};

for k = 1:numel(requiredLinear)
    if ~isprop(out,requiredLinear{k})
        error('Variabile To Workspace lineare mancante in out: %s', ...
              requiredLinear{k});
    end
end

[tAlpha,alpha] = localTimeAndData(out.(names.alpha));
[tBeta,beta]   = localTimeAndData(out.(names.beta));
[tUcmd,uCmd]   = localTimeAndData(out.(names.uCmd));
[tUsat,uSat]   = localTimeAndData(out.(names.uSat));

%% ------------------------------------------------------------
% NON LINEARE
% ------------------------------------------------------------

requiredNL = { ...
    names.alphaNL, ...
    names.betaNL, ...
    names.uCmdNL, ...
    names.uSatNL};

for k = 1:numel(requiredNL)
    if ~isprop(out,requiredNL{k})
        error('Variabile To Workspace non lineare mancante in out: %s', ...
              requiredNL{k});
    end
end

[tAlphaNL,alphaNL] = localTimeAndData(out.(names.alphaNL));
[tBetaNL,betaNL]   = localTimeAndData(out.(names.betaNL));
[tUcmdNL,uCmdNL]   = localTimeAndData(out.(names.uCmdNL));
[tUsatNL,uSatNL]   = localTimeAndData(out.(names.uSatNL));

%% ------------------------------------------------------------
% MISURE LINEARI
% ------------------------------------------------------------

hasMeasured = isprop(out,names.alphaMeas) && ...
              isprop(out,names.betaMeas);

if hasMeasured
    [tAlphaMeas,alphaMeas] = ...
        localTimeAndData(out.(names.alphaMeas));
    [tBetaMeas,betaMeas] = ...
        localTimeAndData(out.(names.betaMeas));
else
    tAlphaMeas = [];
    alphaMeas = [];
    tBetaMeas = [];
    betaMeas = [];
end

%% ------------------------------------------------------------
% MISURE NON LINEARI
% ------------------------------------------------------------

hasMeasuredNL = isprop(out,names.alphaMeasNL) && ...
                isprop(out,names.betaMeasNL);

if hasMeasuredNL
    [tAlphaMeasNL,alphaMeasNL] = ...
        localTimeAndData(out.(names.alphaMeasNL));
    [tBetaMeasNL,betaMeasNL] = ...
        localTimeAndData(out.(names.betaMeasNL));
else
    tAlphaMeasNL = [];
    alphaMeasNL = [];
    tBetaMeasNL = [];
    betaMeasNL = [];
end

%% 5. RIFERIMENTI
% Se il profilo ref e' nel Base Workspace, usiamo quello del modello.
% Altrimenti fallback coerente con il test standard: 0 -> 3 deg su alpha,
% beta = 0. Il tempo del gradino viene preso da ref.alpha.time se presente.
[alphaRef,betaRef,refStepTime] = localReference(tAlpha,tBeta);

%% 6. ERRORI

% Errore rispetto al riferimento
alphaErr = alphaRef - alpha;
betaErr  = betaRef - beta;

% Errore di linearizzazione: NON LINEARE - LINEARE
if numel(tAlpha) ~= numel(tAlphaNL) || any(tAlpha ~= tAlphaNL)
    error('Asse temporale alpha lineare/nonlineare non coincidente.');
end

if numel(tBeta) ~= numel(tBetaNL) || any(tBeta ~= tBetaNL)
    error('Asse temporale beta lineare/nonlineare non coincidente.');
end
alphaLinError = alphaNL - alpha;
betaLinError  = betaNL - beta;

%% 7. METRICHE
metrics = struct();
metrics.controller = simController;
metrics.refStepTime = refStepTime;
metrics.peakErrorAlpha = max(abs(alphaErr));
metrics.peakErrorBeta  = max(abs(betaErr));
metrics.uCmdPeak1      = localPeakColumn(uCmd,1);
metrics.uCmdPeak2      = localPeakColumn(uCmd,2);
metrics.uSatPeak1      = localPeakColumn(uSat,1);
metrics.uSatPeak2      = localPeakColumn(uSat,2);

%% 8. SNAPSHOT STANDARDIZZATO

snapshot = struct();

snapshot.controller = simController;

% =========================
% LINEARE
% =========================
snapshot.tAlpha = tAlpha;
snapshot.alpha = alpha;
snapshot.alphaRef = alphaRef;

snapshot.tBeta = tBeta;
snapshot.beta = beta;
snapshot.betaRef = betaRef;

snapshot.tUcmd = tUcmd;
snapshot.u_cmd = uCmd;

snapshot.tUsat = tUsat;
snapshot.u_sat = uSat;

snapshot.alphaErr = alphaErr;
snapshot.betaErr = betaErr;

% =========================
% NON LINEARE
% =========================
snapshot.tAlphaNL = tAlphaNL;
snapshot.alphaNL = alphaNL;

snapshot.tBetaNL = tBetaNL;
snapshot.betaNL = betaNL;

snapshot.tUcmdNL = tUcmdNL;
snapshot.u_cmdNL = uCmdNL;

snapshot.tUsatNL = tUsatNL;
snapshot.u_satNL = uSatNL;

% =========================
% ERRORE LINEARIZZAZIONE
% =========================
snapshot.alphaLinError = alphaLinError;
snapshot.betaLinError = betaLinError;

% =========================
% MISURE LINEARI
% =========================
snapshot.hasMeasured = hasMeasured;

if hasMeasured
    snapshot.tAlphaMeas = tAlphaMeas;
    snapshot.alphaMeas = alphaMeas;

    snapshot.tBetaMeas = tBetaMeas;
    snapshot.betaMeas = betaMeas;
end

% =========================
% MISURE NON LINEARI
% =========================
snapshot.hasMeasuredNL = hasMeasuredNL;

if hasMeasuredNL
    snapshot.tAlphaMeasNL = tAlphaMeasNL;
    snapshot.alphaMeasNL = alphaMeasNL;

    snapshot.tBetaMeasNL = tBetaMeasNL;
    snapshot.betaMeasNL = betaMeasNL;
end

snapshot.refStepTime = refStepTime;
snapshot.metrics = metrics;

snapshotFile = fullfile(resultDir, ...
    sprintf('%s.mat',simController));

save(snapshotFile,'snapshot','-v7.3');

fprintf('\n============================================================\n');
fprintf('RISULTATO SIMULINK SALVATO\n');
fprintf('Controller : %s\n',upper(simController));
fprintf('File       : %s\n',snapshotFile);
fprintf('============================================================\n');
fprintf('Reference step time = %.3f s\n',refStepTime);
fprintf('Peak |e_alpha| = %.6g rad (%.3f deg)\n', ...
    metrics.peakErrorAlpha,rad2deg(metrics.peakErrorAlpha));
fprintf('Peak |e_beta|  = %.6g rad (%.3f deg)\n', ...
    metrics.peakErrorBeta,rad2deg(metrics.peakErrorBeta));
fprintf('Peak u_cmd     = [%.6g  %.6g]\n',metrics.uCmdPeak1,metrics.uCmdPeak2);
fprintf('Peak u_sat     = [%.6g  %.6g]\n',metrics.uSatPeak1,metrics.uSatPeak2);

%% ========================================================================
% FUNZIONI LOCALI
% ========================================================================
function [alphaRef,betaRef,tStep] = localReference(tAlpha,tBeta)
    alpha0 = 0;
    alphaFinal = deg2rad(3);
    beta0 = 0;
    betaFinal = 0;
    tStep = 1.0;

    if evalin('base',"exist('ref','var')")
        ref = evalin('base','ref');
        if isfield(ref,'alpha')
            if isfield(ref.alpha,'initial'), alpha0 = ref.alpha.initial; end
            if isfield(ref.alpha,'final'),   alphaFinal = ref.alpha.final; end
            if isfield(ref.alpha,'time'),    tStep = ref.alpha.time; end
        end
        if isfield(ref,'beta')
            if isfield(ref.beta,'initial'), beta0 = ref.beta.initial; end
            if isfield(ref.beta,'final'),   betaFinal = ref.beta.final; end
        end
    end

    alphaRef = alpha0*ones(size(tAlpha));
    alphaRef(tAlpha >= tStep) = alphaFinal;
    betaRef  = beta0*ones(size(tBeta));
    betaRef(tBeta >= tStep) = betaFinal;
end

function [t,y] = localTimeAndData(sig)
    if isa(sig,'timeseries')
        t = sig.Time(:);
        y = sig.Data;
    elseif isstruct(sig)
        if isfield(sig,'time')
            t = sig.time(:);
        elseif isfield(sig,'Time')
            t = sig.Time(:);
        else
            t = [];
        end
        if isfield(sig,'signals') && isfield(sig.signals,'values')
            y = sig.signals.values;
        elseif isfield(sig,'Signals') && isfield(sig.Signals,'values')
            y = sig.Signals.values;
        elseif isfield(sig,'values')
            y = sig.values;
        elseif isfield(sig,'Values')
            y = sig.Values;
        else
            error('Formato struct del segnale non riconosciuto.');
        end
    elseif isnumeric(sig)
        y = sig;
        t = (0:size(y,1)-1).';
        warning('Segnale numerico senza tempo: uso indice campione come asse temporale.');
    else
        error('Tipo di segnale non supportato: %s',class(sig));
    end

    y = squeeze(y);
    if isempty(t)
        t = (0:size(y,1)-1).';
    end
    t = t(:);

    if isvector(y)
        y = y(:);
    end
    if size(y,1) ~= numel(t) && size(y,2) == numel(t)
        y = y.';
    end
    if size(y,1) ~= numel(t)
        error('Numero campioni tempo/segnale incompatibile: %d vs %d.',numel(t),size(y,1));
    end
end

function p = localPeakColumn(y,col)
    if size(y,2) < col
        p = NaN;
    else
        p = max(abs(y(:,col)));
    end
end
