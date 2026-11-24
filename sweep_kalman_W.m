%% ========================================================================
% SWEEP DEL PESO DI PROCESSO W DEL FILTRO DI KALMAN
% (LQG 1-DOF, LQG 2-DOF, LQGI)
%
% Per ogni combinazione (W_alpha, W_beta) nella griglia Wgrid:
%   1) rigenera i tre controllori (Kalman incluso) con quel W
%      -> LQG_reg_1DOF.m, LQG_reg_2DOF.m, LQG_con_integratore.m
%   2) lancia la simulazione Simulink (lineare + non lineare, tutti e
%      tre i controllori, esattamente come fatto finora a mano)
%   3) estrae SOLO le metriche numeriche (niente grafici) con
%      estrai_metriche_LQG.m
%
% Cosi' non serve piu' rieseguire a mano i 3 controllori + Simulink +
% analisi_LQG.m per ogni valore di W: lo fa questo script.
%
% PRE-REQUISITI:
%   - I 3 script LQG_reg_1DOF.m, LQG_reg_2DOF.m, LQG_con_integratore.m
%     devono avere la guardia "if ~exist('W','var') ... end" al posto
%     della riga "W = diag([0.01 0.01])" originale (gia' applicata nei
%     file allegati), altrimenti sovrascrivono sempre W = 0.01I e lo
%     sweep non ha effetto.
%   - dataset_elicottero.m deve essere gia' stato eseguito UNA VOLTA in
%     questa sessione (stesso base workspace):
%         run('dataset_elicottero.m')
% ========================================================================

clc;

if ~exist('A_ext','var')
    error(['Esegui prima dataset_elicottero.m in questa sessione ' ...
           '(base workspace) prima di lanciare lo sweep.']);
end

modelName = 'sim_elicottero';
tipo_test = 2;   % stessa convenzione usata in analisi_LQG.m

GkChoice = [B_nom; zeros(6,2)];

%% --- Griglia dei pesi da testare: colonne = [W_alpha, W_beta] ----------
% Punto di partenza: qualche valore simmetrico (come gia' provato a
% mano: 0.001, 0.01, ...) piu' alcune combinazioni ASIMMETRICHE, dato
% che alpha e beta hanno un rapporto segnale/rumore molto diverso.
Wgrid = [ ...
    0.001 0.001;
    0.003 0.003;
    0.010 0.010;   
    0.030 0.030;
    0.100 0.100;
    0.010 0.003;
    0.010 0.001;
    0.030 0.010];

nSweep = size(Wgrid,1);

sweepResults = table( ...
    'Size',[nSweep 10], ...
    'VariableTypes',repmat({'double'},1,10), ...
    'VariableNames',{ ...
        'W_alpha','W_beta', ...
        'RMSest_alpha_LQG2','RMSest_beta_LQG2', ...
        'RMSest_alpha_LQGI','RMSest_beta_LQGI', ...
        'Trec_beta_LQG2','Trec_beta_LQGI', ...
        'PeakBetaDist_LQG2','PeakBetaDist_LQGI'});

set(0,'DefaultFigureVisible','off');   % niente figure a schermo durante lo sweep

for k = 1:nSweep

    W = diag(Wgrid(k,:));

    run('LQG_reg_1DOF.m');
    run('LQG_reg_2DOF.m');
    run('LQG_con_integratore.m');

   % ------------------------------------------------------------
    % SIMULAZIONE
    % ------------------------------------------------------------
    
    clear out
    
    simOut = sim(modelName, ...
        'ReturnWorkspaceOutputs','on');
    disp(simOut.who)   % TEMPORANEAMENTE, per capire la struttura
    % ------------------------------------------------------------
    % RECUPERO RISULTATI
    % ------------------------------------------------------------
    
    if isa(simOut,'Simulink.SimulationOutput')
    
        nomi = simOut.who;
    
        % Caso 1:
        % Simulink ha restituito una variabile "out"
        if any(strcmp(nomi,'out'))
    
            out = simOut.get('out');
    
        else
    
            % Caso 2:
            % i segnali sono direttamente dentro SimulationOutput
            out = simOut;
    
        end
    
    else
    
        out = simOut;
    
    end
    
    % ------------------------------------------------------------
    % CONTROLLO CHE I SEGNALI NECESSARI ESISTANO
    % ------------------------------------------------------------
    
    segnaliNecessari = { ...
        'alpha_LQG2', ...
        'beta_LQG2', ...
        'xhat_LQG2', ...
        'alpha_LQGI', ...
        'beta_LQGI', ...
        'xhat_LQGI', ...
        'alpha_LQG2_nl', ...
        'beta_LQG2_nl', ...
        'xhat_LQG2_nl', ...
        'alpha_LQGI_nl', ...
        'beta_LQGI_nl', ...
        'xhat_LQGI_nl'};
    
    mancanti = {};
    
    for ii = 1:numel(segnaliNecessari)
    
        nomeSegnale = segnaliNecessari{ii};
    
        try
            out.(nomeSegnale);
        catch
            mancanti{end+1} = nomeSegnale;
        end
    
    end
    
    if ~isempty(mancanti)
    
        disp('Segnali disponibili in simOut:')
        disp(simOut.who)
    
        disp('Segnali mancanti:')
        disp(mancanti)
    
        error('La simulazione non ha prodotto tutti i segnali necessari.')
    
    end

    m = estrai_metriche_LQG(out, tipo_test);

    sweepResults.W_alpha(k) = Wgrid(k,1);
    sweepResults.W_beta(k)  = Wgrid(k,2);

    sweepResults.RMSest_alpha_LQG2(k) = m.LQG2.RMS_est_alpha_lin;
    sweepResults.RMSest_beta_LQG2(k)  = m.LQG2.RMS_est_beta_lin;
    sweepResults.RMSest_alpha_LQGI(k) = m.LQGI.RMS_est_alpha_lin;
    sweepResults.RMSest_beta_LQGI(k)  = m.LQGI.RMS_est_beta_lin;

    sweepResults.Trec_beta_LQG2(k) = m.LQG2.ris_lin.T_rec_beta;
    sweepResults.Trec_beta_LQGI(k) = m.LQGI.ris_lin.T_rec_beta;

    sweepResults.PeakBetaDist_LQG2(k) = m.LQG2.ris_lin.peak_beta_dist;
    sweepResults.PeakBetaDist_LQGI(k) = m.LQGI.ris_lin.peak_beta_dist;

    close all;  % chiude anche le figure create da calcola_prestazioni_LQG

    fprintf(['[%d/%d] W_alpha=%.4f  W_beta=%.4f  ->  ' ...
             'RMSest_beta(LQGI)=%.5f rad   Trec_beta(LQGI)=%.3f s\n'], ...
        k, nSweep, Wgrid(k,1), Wgrid(k,2), ...
        sweepResults.RMSest_beta_LQGI(k), sweepResults.Trec_beta_LQGI(k));

end



save('sweep_kalman_results.mat','sweepResults','Wgrid');
disp(sweepResults);

%% --- Vista rapida del trade-off rumore/reiezione (canale beta, LQGI) --
figure('Name','Pareto Kalman W - LQGI (beta)');
scatter(sweepResults.RMSest_beta_LQGI, sweepResults.Trec_beta_LQGI, 60, 'filled');
for k = 1:nSweep
    text(sweepResults.RMSest_beta_LQGI(k), sweepResults.Trec_beta_LQGI(k), ...
        sprintf('  (%.3f, %.3f)', Wgrid(k,1), Wgrid(k,2)), 'FontSize', 8);
end
grid on
xlabel('RMS errore di stima \beta [rad]  (proxy filtraggio rumore, piu'' basso = meglio)')
ylabel('Tempo di recupero \beta post-disturbo [s]  (proxy reiezione, piu'' basso = meglio)')
title('Trade-off Kalman W - LQGI - ogni punto = (W_\alpha, W_\beta)')
