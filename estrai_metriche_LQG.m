function metriche = estrai_metriche_LQG(out, tipo_test)
%ESTRAI_METRICHE_LQG Estrae le metriche numeriche di tracking e di stima
%   Kalman dai segnali contenuti in "out" (l'output di sim_elicottero),
%   SENZA generare i grafici di analisi_LQG.m.
%
%   Pensata per essere richiamata dentro un ciclo di sweep parametrico
%   (vedi sweep_kalman_W.m): analisi_LQG.m e' troppo pesante da rilanciare
%   ad ogni iterazione perche' genera decine di figure.
%
%   metriche = estrai_metriche_LQG(out, tipo_test)
%
%   tipo_test : 1, 2 o 3, stessa convenzione di analisi_LQG.m
%               (default: 2)
%
%   Riusa calcola_prestazioni_LQG.m cosi' com'e' (stesse metriche di
%   tracking/disturbance rejection gia' usate nel resto della tesi),
%   quindi deve trovarsi sul path.
%
%   Restituisce una struct con un campo per ciascun controllore
%   analizzato ('LQG','LQG2' e/o 'LQGI', a seconda di tipo_test), ognuno
%   con:
%       .ris_lin, .ris_nl          -> output di calcola_prestazioni_LQG
%                                      (modello lineare / non lineare)
%       .RMS_est_alpha_lin/beta    -> errore di stima Kalman (x - xhat),
%                                      modello lineare
%       .RMS_est_alpha_nl/beta     -> idem, modello non lineare

if nargin < 2
    tipo_test = 2;
end

controllers     = {'LQG','LQG2','LQGI'};
controllerNames = {'LQG 1-DOF','LQG 2-DOF','LQGI'};

% Stessa convenzione di analisi_LQG.m
if tipo_test == 1 || tipo_test == 2
    idx_confronto = [2 3];   % LQG 2-DOF, LQGI
else
    idx_confronto = [1 2 3]; % LQG 1-DOF, LQG 2-DOF, LQGI
end

metriche = struct();

for j = 1:numel(idx_confronto)
    i    = idx_confronto(j);
    suff = controllers{i};

    %% --- Segnali dal workspace di simulazione ---------------------
    alpha   = out.(['alpha_' suff]);
    beta    = out.(['beta_'  suff]);
    xhat    = out.(['xhat_'  suff]);

    alphaNL = out.(['alpha_' suff '_nl']);
    betaNL  = out.(['beta_'  suff '_nl']);
    xhatNL  = out.(['xhat_'  suff '_nl']);

    ris_lin = calcola_prestazioni_LQG( ...
        alpha, ...
        beta, ...
        tipo_test, ...
        [controllerNames{i} ' - Lineare'], ...
        false);
    
    ris_nl = calcola_prestazioni_LQG( ...
        alphaNL, ...
        betaNL, ...
        tipo_test, ...
        [controllerNames{i} ' - Non lineare'], ...
        false);

    %% --- Errore di stima Kalman: MODELLO LINEARE -------------------
    tAlpha = alpha.Time(:);  alphaV = alpha.Data(:);
    tBeta  = beta.Time(:);   betaV  = beta.Data(:);

    tXhat  = xhat.Time(:);
    xhatV  = squeeze(xhat.Data);
    if size(xhatV,2) ~= 8 && size(xhatV,1) == 8
        xhatV = xhatV.';
    end
    if size(xhatV,2) ~= 8
        error('xhat_%s deve contenere 8 colonne.', suff)
    end

    alpha_real_xhat = interp1(tAlpha, alphaV, tXhat, 'linear', 'extrap');
    beta_real_xhat  = interp1(tBeta,  betaV,  tXhat, 'linear', 'extrap');

    errEstAlpha = alpha_real_xhat - xhatV(:,1);
    errEstBeta  = beta_real_xhat  - xhatV(:,3);

    %% --- Errore di stima Kalman: MODELLO NON LINEARE ----------------
    tAlphaNL = alphaNL.Time(:);  alphaNLV = alphaNL.Data(:);
    tBetaNL  = betaNL.Time(:);   betaNLV  = betaNL.Data(:);

    tXhatNL  = xhatNL.Time(:);
    xhatNLV  = squeeze(xhatNL.Data);
    if size(xhatNLV,2) ~= 8 && size(xhatNLV,1) == 8
        xhatNLV = xhatNLV.';
    end
    if size(xhatNLV,2) ~= 8
        error('xhat_%s_nl deve contenere 8 colonne.', suff)
    end

    alpha_real_xhatNL = interp1(tAlphaNL, alphaNLV, tXhatNL, 'linear', 'extrap');
    beta_real_xhatNL  = interp1(tBetaNL,  betaNLV,  tXhatNL, 'linear', 'extrap');

    errEstAlphaNL = alpha_real_xhatNL - xhatNLV(:,1);
    errEstBetaNL  = beta_real_xhatNL  - xhatNLV(:,3);

    %% --- Raccolta risultati -----------------------------------------
    m = struct();
    m.ris_lin = ris_lin;
    m.ris_nl  = ris_nl;

    m.RMS_est_alpha_lin = rms(errEstAlpha);
    m.RMS_est_beta_lin  = rms(errEstBeta);
    m.RMS_est_alpha_nl  = rms(errEstAlphaNL);
    m.RMS_est_beta_nl   = rms(errEstBetaNL);

    metriche.(suff) = m;
end

end
