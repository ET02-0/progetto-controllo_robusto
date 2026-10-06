close all
clc

% ================================================================
% IMPOSTAZIONE GRAFICA PER LA RELAZIONE
% ================================================================

set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k', ...
    'defaultAxesYColor','k', ...
    'defaultAxesZColor','k', ...
    'defaultAxesGridColor',[0.75 0.75 0.75], ...
    'defaultTextColor','k', ...
    'defaultAxesBox','on', ...
    'defaultAxesFontSize',11);

%% ================================================================
% ANALISI PRESTAZIONI LQG / LQGI - SIMULINK
%
% I controllori analizzati sono sempre:
%
%   1) LQG 1-DOF
%   2) LQG 2-DOF senza integratore
%   3) LQG 2-DOF con integratore (LQGI)
%
% Tutti i segnali lineari e non lineari sono già presenti
% contemporaneamente in "out".
%
% CONVENZIONE DEI NOMI:
%
% LQG 1-DOF:
%   alpha_LQG
%   beta_LQG
%   alpha_meas_LQG
%   beta_meas_LQG
%   u_cmd_LQG
%   u_sat_LQG
%   delta_F_LQG
%
% LQG 1-DOF non lineare:
%   alpha_LQG_nl
%   beta_LQG_nl
%   alpha_meas_LQG_nl
%   beta_meas_LQG_nl
%   u_cmd_LQG_nl
%   u_sat_LQG_nl
%   delta_F_LQG_nl
%
% LQG 2-DOF:
%   ..._LQG2
%   ..._LQG2_nl
%
% LQGI:
%   ..._LQGI
%   ..._LQGI_nl
%
% tipo_test:
%   1 = tracking nominale, senza rumore/disturbo
%   2 = tracking + disturbo + rumore
%   3 = disturbo + rumore con riferimento nullo
%% ================================================================

tipo_test = 2;

%% ================================================================
% CONTROLLORI DA CONFRONTARE NEL TEST CORRENTE
%
% Tipo 1: tracking nominale -> LQG 2-DOF + LQGI
% Tipo 2: tracking + disturbo + rumore -> LQG 2-DOF + LQGI
% Tipo 3: disturbo + rumore, riferimento nullo
%          -> LQG 1-DOF + LQG 2-DOF + LQGI
%% ================================================================

idxLQG1 = 1;
idxLQG2 = 2;
idxLQGI = 3;

if tipo_test == 1 || tipo_test == 2
    idx_confronto = [idxLQG2 idxLQGI];
else
    idx_confronto = [idxLQG1 idxLQG2 idxLQGI];
end


%% ================================================================
% 1. RIFERIMENTI
%% ================================================================

switch tipo_test

    case 1
        r_alpha = deg2rad(3);
        r_beta  = 0;

    case 2
        r_alpha = deg2rad(3);
        r_beta  = 0;

    case 3
        r_alpha = 0;
        r_beta  = 0;

    otherwise
        error('tipo_test deve essere 1, 2 oppure 3.')

end

%% ================================================================
% 2. PUNTO DI EQUILIBRIO E RIFERIMENTI ASSOLUTI
%
% I segnali alpha e beta presenti in "out" sono espressi come
% deviazioni dall'equilibrio:
%
%   alpha = Delta_alpha
%   beta  = Delta_beta
%
% Per i grafici delle grandezze assolute si ricostruiscono quindi:
%
%   alpha_abs = alpha_0 + Delta_alpha
%   beta_abs  = beta_0  + Delta_beta
%% ================================================================

if ~exist('alpha_0','var') || ~exist('beta_0','var')
    error(['Le variabili alpha_0 e beta_0 non sono presenti nel ', ...
           'workspace. Definirle prima di eseguire analisi_LQG.'])
end

r_alpha_abs = alpha_0 + r_alpha;
r_beta_abs  = beta_0  + r_beta;


%% ================================================================
% 2. DEFINIZIONE DEI TRE CONTROLLORI
%% ================================================================

controllers = { ...
    'LQG', ...
    'LQG2', ...
    'LQGI'};

controllerNames = { ...
    'LQG 1-DOF', ...
    'LQG 2-DOF', ...
    'LQGI'};

Nc = numel(controllers);
labelsConfronto = controllerNames(idx_confronto);


%% ================================================================
% 3. LETTURA DI TUTTI I SEGNALI DA "out"
%% ================================================================
D = struct();
for i = 1:Nc

    suff = controllers{i};

    %% ------------------------------------------------------------
    % MODELLO LINEARE
    %% ------------------------------------------------------------

    D(i).alpha        = out.(['alpha_' suff]);
    D(i).beta         = out.(['beta_' suff]);

    D(i).alpha_meas   = out.(['alpha_meas_' suff]);
    D(i).beta_meas    = out.(['beta_meas_' suff]);

    D(i).u_cmd        = out.(['u_cmd_' suff]);
    D(i).u_sat        = out.(['u_sat_' suff]);

    D(i).delta_F      = out.(['delta_F_' suff]);

    %% ------------------------------------------------------------
    % STIMA DELLO STATO - FILTRO DI KALMAN
    %
    % xhat contiene:
    % [alpha_hat, alphadot_hat, beta_hat, betadot_hat,
    %  act1_hat, act1dot_hat, act2_hat, act2dot_hat]
    %% ------------------------------------------------------------

    D(i).xhat = out.(['xhat_' suff]);


    %% ------------------------------------------------------------
    % MODELLO NON LINEARE
    %% ------------------------------------------------------------

    D(i).alphaNL      = out.(['alpha_' suff '_nl']);
    D(i).betaNL       = out.(['beta_' suff '_nl']);

    D(i).alpha_measNL = out.(['alpha_meas_' suff '_nl']);
    D(i).beta_measNL  = out.(['beta_meas_' suff '_nl']);

    D(i).u_cmdNL      = out.(['u_cmd_' suff '_nl']);
    D(i).u_satNL      = out.(['u_sat_' suff '_nl']);

    D(i).delta_FNL    = out.(['delta_F_' suff '_nl']);

    D(i).xhatNL = out.(['xhat_' suff '_nl']);


    %% ------------------------------------------------------------
    % VETTORI NUMERICI LINEARI
    %% ------------------------------------------------------------

    D(i).tAlpha = D(i).alpha.Time(:);
    D(i).alphaV = D(i).alpha.Data(:);

    D(i).tBeta = D(i).beta.Time(:);
    D(i).betaV = D(i).beta.Data(:);

    % Grandezze assolute per la visualizzazione
    D(i).alphaAbsV = alpha_0 + D(i).alphaV;
    D(i).betaAbsV  = beta_0  + D(i).betaV;

    D(i).tAlphaMeas = D(i).alpha_meas.Time(:);
    D(i).alphaMeasV = D(i).alpha_meas.Data(:);

    D(i).tBetaMeas = D(i).beta_meas.Time(:);
    D(i).betaMeasV = D(i).beta_meas.Data(:);

    % Misure assolute per la visualizzazione
    D(i).alphaMeasAbsV = alpha_0 + D(i).alphaMeasV;
    D(i).betaMeasAbsV  = beta_0  + D(i).betaMeasV;

    D(i).tUcmd = D(i).u_cmd.Time(:);
    D(i).uCmdV = D(i).u_cmd.Data;

    D(i).tUsat = D(i).u_sat.Time(:);
    D(i).uSatV = D(i).u_sat.Data;

    D(i).tDeltaF = D(i).delta_F.Time(:);
    D(i).deltaFV = D(i).delta_F.Data;

    %% ------------------------------------------------------------
    % STIMA KALMAN - MODELLO LINEARE
    %% ------------------------------------------------------------

    D(i).tXhat = D(i).xhat.Time(:);
    D(i).xhatV = squeeze(D(i).xhat.Data);

    % Garantisce la convenzione N x 8
    if size(D(i).xhatV,2) ~= 8 && size(D(i).xhatV,1) == 8
        D(i).xhatV = D(i).xhatV.';
    end

    if size(D(i).xhatV,2) ~= 8
        error('xhat_%s deve contenere 8 colonne.', suff)
    end


    %% ------------------------------------------------------------
    % VETTORI NUMERICI NON LINEARI
    %% ------------------------------------------------------------

    D(i).tAlphaNL = D(i).alphaNL.Time(:);
    D(i).alphaNLV = D(i).alphaNL.Data(:);

    D(i).tBetaNL = D(i).betaNL.Time(:);
    D(i).betaNLV = D(i).betaNL.Data(:);
    % Grandezze assolute per la visualizzazione
    D(i).alphaNLAbsV = alpha_0 + D(i).alphaNLV;
    D(i).betaNLAbsV  = beta_0  + D(i).betaNLV;

    D(i).tAlphaMeasNL = D(i).alpha_measNL.Time(:);
    D(i).alphaMeasNLV = D(i).alpha_measNL.Data(:);

    D(i).tBetaMeasNL = D(i).beta_measNL.Time(:);
    D(i).betaMeasNLV = D(i).beta_measNL.Data(:);

    % Misure assolute per la visualizzazione
    D(i).alphaMeasNLAbsV = alpha_0 + D(i).alphaMeasNLV;
    D(i).betaMeasNLAbsV  = beta_0  + D(i).betaMeasNLV;

    D(i).tUcmdNL = D(i).u_cmdNL.Time(:);
    D(i).uCmdNLV = D(i).u_cmdNL.Data;

    D(i).tUsatNL = D(i).u_satNL.Time(:);
    D(i).uSatNLV = D(i).u_satNL.Data;

    D(i).tDeltaFNL = D(i).delta_FNL.Time(:);
    D(i).deltaFNLV = D(i).delta_FNL.Data;

    %% ------------------------------------------------------------
    % STIMA KALMAN - MODELLO NON LINEARE
    %% ------------------------------------------------------------

    D(i).tXhatNL = D(i).xhatNL.Time(:);
    D(i).xhatNLV = squeeze(D(i).xhatNL.Data);

    % Garantisce la convenzione N x 8
    if size(D(i).xhatNLV,2) ~= 8 && size(D(i).xhatNLV,1) == 8
        D(i).xhatNLV = D(i).xhatNLV.';
    end

    if size(D(i).xhatNLV,2) ~= 8
        error('xhat_%s_nl deve contenere 8 colonne.', suff)
    end


    %% ------------------------------------------------------------
    % CONTROLLO DIMENSIONI
    %% ------------------------------------------------------------

    if size(D(i).uCmdV,2) < 2
        error('u_cmd_%s deve contenere due colonne.', suff)
    end

    if size(D(i).uSatV,2) < 2
        error('u_sat_%s deve contenere due colonne.', suff)
    end

    if size(D(i).deltaFV,2) < 2
        error('delta_F_%s deve contenere due colonne.', suff)
    end

    if size(D(i).uCmdNLV,2) < 2
        error('u_cmd_%s_nl deve contenere due colonne.', suff)
    end

    if size(D(i).uSatNLV,2) < 2
        error('u_sat_%s_nl deve contenere due colonne.', suff)
    end

    if size(D(i).deltaFNLV,2) < 2
        error('delta_F_%s_nl deve contenere due colonne.', suff)
    end

end


%% ================================================================
% 4. DISTURBI AERODINAMICI
%    Ricostruiti direttamente dalla variabile "aero" del dataset
%% ================================================================
if tipo_test == 2 || tipo_test == 3
    hasAero = false;
    
    if evalin('base',"exist('aero','var')")
    
        aero = evalin('base','aero');
    
        if isfield(aero,'enable') && aero.enable ~= 0
    
            % Usiamo il tempo del primo controllore come riferimento
            t_aero = D(idx_confronto(1)).tAlpha;
    
            % Inizializzazione
            dist_aero_alpha = zeros(size(t_aero));
            dist_aero_beta  = zeros(size(t_aero));
    
            % -----------------------------
            % Disturbo su alpha
            % -----------------------------
            if isfield(aero,'alpha') && ...
               isfield(aero.alpha,'time') && ...
               isfield(aero.alpha,'amplitude')
    
                tDa = aero.alpha.time;
                aDa = aero.alpha.amplitude;
    
                dist_aero_alpha(t_aero >= tDa) = aDa;
    
            end
    
            % -----------------------------
            % Disturbo su beta
            % -----------------------------
            if isfield(aero,'beta') && ...
               isfield(aero.beta,'time') && ...
               isfield(aero.beta,'amplitude')
    
                tDb = aero.beta.time;
                aDb = aero.beta.amplitude;
    
                dist_aero_beta(t_aero >= tDb) = aDb;
    
            end
    
            hasAero = true;
    
        end
    
    else
    
        warning('Variabile "aero" non trovata nel Base Workspace.');
    
        t_aero = [];
        dist_aero_alpha = [];
        dist_aero_beta = [];
    
    end
end

%% ================================================================
% 5. ANALISI PRESTAZIONI DEI TRE CONTROLLORI
%% ================================================================

fprintf('\n\n')
fprintf('############################################################\n')
fprintf('# ANALISI PRESTAZIONI LQG / LQGI\n')
fprintf('############################################################\n')

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    fprintf('\n\n')
    fprintf('============================================================\n')
    fprintf(' %s - MODELLO LINEARE\n',controllerNames{i})
    fprintf('============================================================\n')

    D(i).ris_lin = calcola_prestazioni_LQG( ...
        D(i).alpha, ...
        D(i).beta, ...
        tipo_test, ...
        [controllerNames{i} ' - Lineare']);


    fprintf('\n\n')
    fprintf('============================================================\n')
    fprintf(' %s - MODELLO NON LINEARE\n',controllerNames{i})
    fprintf('============================================================\n')

    D(i).ris_nl = calcola_prestazioni_LQG( ...
        D(i).alphaNL, ...
        D(i).betaNL, ...
        tipo_test, ...
        [controllerNames{i} ' - Non lineare']);

end


%% ================================================================
% 6. ERRORE DI LINEARIZZAZIONE
%% ================================================================

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    % Interpolazione del modello lineare sulla griglia non lineare
    D(i).alphaLinInterp = interp1( ...
        D(i).tAlpha, ...
        D(i).alphaV, ...
        D(i).tAlphaNL, ...
        'linear', ...
        'extrap');

    D(i).betaLinInterp = interp1( ...
        D(i).tBeta, ...
        D(i).betaV, ...
        D(i).tBetaNL, ...
        'linear', ...
        'extrap');


    D(i).errAlpha = ...
        D(i).alphaNLV - D(i).alphaLinInterp;

    D(i).errBeta = ...
        D(i).betaNLV - D(i).betaLinInterp;

end

%% ================================================================
% 6B. ERRORE DI STIMA DEL FILTRO DI KALMAN
%
% e_hat = x - x_hat
%
% Si considerano:
%   alpha_hat = xhat(:,1)
%   beta_hat  = xhat(:,3)
%
% I segnali reali vengono interpolati sulla griglia temporale
% dello stimatore per garantire il confronto campione per campione.
%% ================================================================

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    %% ------------------------------------------------------------
    % MODELLO LINEARE
    %% ------------------------------------------------------------

    alpha_real_xhat = interp1( ...
        D(i).tAlpha, ...
        D(i).alphaV, ...
        D(i).tXhat, ...
        'linear', ...
        'extrap');

    beta_real_xhat = interp1( ...
        D(i).tBeta, ...
        D(i).betaV, ...
        D(i).tXhat, ...
        'linear', ...
        'extrap');

    % Stati stimati
    D(i).alphaHatV = D(i).xhatV(:,1);
    D(i).betaHatV  = D(i).xhatV(:,3);

    % Errore di stima: x - xhat
    D(i).errEstAlpha = ...
        alpha_real_xhat - D(i).alphaHatV;

    D(i).errEstBeta = ...
        beta_real_xhat - D(i).betaHatV;


    %% ------------------------------------------------------------
    % MODELLO NON LINEARE
    %% ------------------------------------------------------------

    alpha_real_xhatNL = interp1( ...
        D(i).tAlphaNL, ...
        D(i).alphaNLV, ...
        D(i).tXhatNL, ...
        'linear', ...
        'extrap');

    beta_real_xhatNL = interp1( ...
        D(i).tBetaNL, ...
        D(i).betaNLV, ...
        D(i).tXhatNL, ...
        'linear', ...
        'extrap');

    % Stati stimati
    D(i).alphaHatNLV = D(i).xhatNLV(:,1);
    D(i).betaHatNLV  = D(i).xhatNLV(:,3);

    D(i).alphaHatAbsV   = alpha_0 + D(i).alphaHatV;
    D(i).betaHatAbsV    = beta_0  + D(i).betaHatV;
    
    D(i).alphaHatNLAbsV = alpha_0 + D(i).alphaHatNLV;
    D(i).betaHatNLAbsV  = beta_0  + D(i).betaHatNLV;

    % Errore di stima
    D(i).errEstAlphaNL = ...
        alpha_real_xhatNL - D(i).alphaHatNLV;

    D(i).errEstBetaNL = ...
        beta_real_xhatNL - D(i).betaHatNLV;

end




%% ================================================================
% 7. GRAFICI: LINEARE vs NON LINEARE
%    UNO PER CIASCUN CONTROLLore
%% ================================================================

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    figure( ...
        'Name',[controllerNames{i} ' - Linear vs Nonlinear'], ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(i).tAlpha, ...
        rad2deg(D(i).alphaAbsV), ...
        '--', ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(i).tAlphaNL, ...
        rad2deg(D(i).alphaNLAbsV), ...
        '-', ...
        'LineWidth',1.5)

    yline( ...
        rad2deg(r_alpha_abs), ...
        ':', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\alpha [deg]')

    title([controllerNames{i} ' - Pitch'])

    legend( ...
        'Linearized', ...
        'Nonlinear', ...
        'Reference', ...
        'Location','best')


    subplot(2,1,2)

    plot( ...
        D(i).tBeta, ...
        rad2deg(D(i).betaAbsV), ...
        '--', ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(i).tBetaNL, ...
        rad2deg(D(i).betaNLAbsV), ...
        '-', ...
        'LineWidth',1.5)

    yline( ...
        rad2deg(r_beta_abs), ...
        ':', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\beta [deg]')

    title([controllerNames{i} ' - Yaw'])

    legend( ...
        'Linearized', ...
        'Nonlinear', ...
        'Reference', ...
        'Location','best')

end


%% ================================================================
% 8. GRAFICI: ERRORE DI LINEARIZZAZIONE
%% ================================================================

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    figure( ...
        'Name',[controllerNames{i} ' - Linearization Error'], ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(i).tAlphaNL, ...
        rad2deg(D(i).errAlpha), ...
        'LineWidth',1.5)

    hold on

    yline(0,'k--','LineWidth',0.7)

    grid on

    xlabel('Time [s]')
    ylabel('\alpha_{NL}-\alpha_{LIN} [deg]')

    title([controllerNames{i} ' - Pitch linearization error'])


    subplot(2,1,2)

    plot( ...
        D(i).tBetaNL, ...
        rad2deg(D(i).errBeta), ...
        'LineWidth',1.5)

    hold on

    yline(0,'k--','LineWidth',0.7)

    grid on

    xlabel('Time [s]')
    ylabel('\beta_{NL}-\beta_{LIN} [deg]')

    title([controllerNames{i} ' - Yaw linearization error'])

end


%% ================================================================
% ERRORE DI LINEARIZZAZIONE - CONFRONTO TRA TUTTI I CONTROLLORI
%
% Figura comparativa utilizzata nella relazione.
%% ================================================================


figure( ...
    'Name',sprintf('LQG comparison - Linearization Error - Tipo %d',tipo_test), ...
    'Color','w')

% ------------------------------------------------------------
% PITCH
% ------------------------------------------------------------

subplot(2,1,1)
hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tAlphaNL, ...
        rad2deg(D(i).errAlpha), ...
        'LineWidth',1.5, ...
        'DisplayName',controllerNames{i})

end

yline(0,'k--','LineWidth',0.8)

grid on

xlabel('Time [s]')
ylabel('\alpha_{NL}-\alpha_{LIN} [deg]')

title('Errore di linearizzazione - Pitch')

legend('Location','best')


% ------------------------------------------------------------
% YAW
% ------------------------------------------------------------

subplot(2,1,2)
hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tBetaNL, ...
        rad2deg(D(i).errBeta), ...
        'LineWidth',1.5, ...
        'DisplayName',controllerNames{i})

end

yline(0,'k--','LineWidth',0.8)

grid on

xlabel('Time [s]')
ylabel('\beta_{NL}-\beta_{LIN} [deg]')

title('Errore di linearizzazione - Yaw')

legend('Location','best')



%% ================================================================
% 9. COMANDO DI CONTROLLO + SATURAZIONE + USCITA ATTUATORE
%
% Per il modello lineare.
%% ================================================================

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    figure( ...
        'Name',[controllerNames{i} ' - Control Signals'], ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(i).tUcmd, ...
        D(i).uCmdV(:,1), ...
        '--', ...
        'LineWidth',1.3)

    hold on

    plot( ...
        D(i).tUsat, ...
        D(i).uSatV(:,1), ...
        '-', ...
        'LineWidth',1.3)

    plot( ...
        D(i).tDeltaF, ...
        D(i).deltaFV(:,1), ...
        '-.', ...
        'LineWidth',1.3)

    grid on

    xlabel('Time [s]')
    ylabel('\Delta F_1 [N]')

    title([controllerNames{i} ' - Main rotor'])

    legend( ...
        'Command', ...
        'Saturated', ...
        'Actual actuator output', ...
        'Location','best')


    subplot(2,1,2)

    plot( ...
        D(i).tUcmd, ...
        D(i).uCmdV(:,2), ...
        '--', ...
        'LineWidth',1.3)

    hold on

    plot( ...
        D(i).tUsat, ...
        D(i).uSatV(:,2), ...
        '-', ...
        'LineWidth',1.3)

    plot( ...
        D(i).tDeltaF, ...
        D(i).deltaFV(:,2), ...
        '-.', ...
        'LineWidth',1.3)

    grid on

    xlabel('Time [s]')
    ylabel('\Delta F_2 [N]')

    title([controllerNames{i} ' - Tail rotor'])

    legend( ...
        'Command', ...
        'Saturated', ...
        'Actual actuator output', ...
        'Location','best')

end


%% ================================================================
% 10. SEGNALI DI CONTROLLO NON LINEARI
%% ================================================================

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    figure( ...
        'Name',[controllerNames{i} ' - Nonlinear Control Signals'], ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(i).tUcmdNL, ...
        D(i).uCmdNLV(:,1), ...
        '--', ...
        'LineWidth',1.3)

    hold on

    plot( ...
        D(i).tUsatNL, ...
        D(i).uSatNLV(:,1), ...
        '-', ...
        'LineWidth',1.3)

    plot( ...
        D(i).tDeltaFNL, ...
        D(i).deltaFNLV(:,1), ...
        '-.', ...
        'LineWidth',1.3)

    grid on

    xlabel('Time [s]')
    ylabel('\Delta F_1 [N]')

    title([controllerNames{i} ' - Nonlinear main rotor'])

    legend( ...
        'Command', ...
        'Saturated', ...
        'Actual actuator output', ...
        'Location','best')


    subplot(2,1,2)

    plot( ...
        D(i).tUcmdNL, ...
        D(i).uCmdNLV(:,2), ...
        '--', ...
        'LineWidth',1.3)

    hold on

    plot( ...
        D(i).tUsatNL, ...
        D(i).uSatNLV(:,2), ...
        '-', ...
        'LineWidth',1.3)

    plot( ...
        D(i).tDeltaFNL, ...
        D(i).deltaFNLV(:,2), ...
        '-.', ...
        'LineWidth',1.3)

    grid on

    xlabel('Time [s]')
    ylabel('\Delta F_2 [N]')

    title([controllerNames{i} ' - Nonlinear tail rotor'])

    legend( ...
        'Command', ...
        'Saturated', ...
        'Actual actuator output', ...
        'Location','best')

end


%% ================================================================
% 11. MISURATO vs REALE
%    Il segnale misurato è rumoroso: lo visualizziamo come punti
%    e disegniamo il segnale reale sopra, ben visibile.
%% ================================================================
if tipo_test == 2 || tipo_test == 3
    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);
    
        figure( ...
            'Name',[controllerNames{i} ' - Measured vs Real'], ...
            'Color','w')
    
    
        %% ------------------------------------------------------------
        % PITCH
        %% ------------------------------------------------------------
    
        subplot(2,1,1)
    
        % Misurato: punti piccoli, così il rumore non copre il reale
        plot( ...
            D(i).tAlphaMeas, ...
            rad2deg(D(i).alphaMeasAbsV), ...
            '.', ...
            'MarkerSize',3)
    
        hold on
    
        % Reale: linea continua disegnata per ultima
        plot( ...
            D(i).tAlpha, ...
            rad2deg(D(i).alphaAbsV), ...
            'k-', ...
            'LineWidth',2)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('\alpha [deg]')
    
        title([controllerNames{i} ' - Pitch: real vs measured'])
    
        legend( ...
            'Measured', ...
            'Real', ...
            'Location','best')
    
    
        %% ------------------------------------------------------------
        % YAW
        %% ------------------------------------------------------------
    
        subplot(2,1,2)
    
        % Misurato: punti piccoli
        plot( ...
            D(i).tBetaMeas, ...
            rad2deg(D(i).betaMeasAbsV), ...
            '.', ...
            'MarkerSize',3)
    
        hold on
    
        % Reale: linea continua sopra il rumore
        plot( ...
            D(i).tBeta, ...
            rad2deg(D(i).betaAbsV), ...
            'k-', ...
            'LineWidth',2)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('\beta [deg]')
    
        title([controllerNames{i} ' - Yaw: real vs measured'])
    
        legend( ...
            'Measured', ...
            'Real', ...
            'Location','best')
    
    end
end
%% ================================================================
% 11B. ERRORE DI MISURA - UNO PER CIASCUN CONTROLLore
%
% e_alpha = alpha_misurata - alpha_reale
% e_beta  = beta_misurata  - beta_reale
%% ================================================================
if tipo_test == 2 || tipo_test == 3
    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);
    
        % -------------------------------------------------------------
        % Errore di misura
        % -------------------------------------------------------------
        eAlphaMeas = D(i).alphaMeasV - D(i).alphaV;
        eBetaMeas  = D(i).betaMeasV  - D(i).betaV;
    
    
        figure( ...
            'Name',[controllerNames{i} ' - Measurement Error'], ...
            'Color','w')
    
    
        % -------------------------------------------------------------
        % PITCH
        % -------------------------------------------------------------
        subplot(2,1,1)
    
        plot( ...
            D(i).tAlpha, ...
            rad2deg(eAlphaMeas), ...
            'LineWidth',1.3)
    
        hold on
    
        yline(0,'k--','LineWidth',0.8)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('e_\alpha [deg]')
    
        title([controllerNames{i} ' - Pitch measurement error'])
    
        legend( ...
            'Measured - Real', ...
            'Zero error', ...
            'Location','best')
    
    
        % -------------------------------------------------------------
        % YAW
        % -------------------------------------------------------------
        subplot(2,1,2)
    
        plot( ...
            D(i).tBeta, ...
            rad2deg(eBetaMeas), ...
            'LineWidth',1.3)
    
        hold on
    
        yline(0,'k--','LineWidth',0.8)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('e_\beta [deg]')
    
        title([controllerNames{i} ' - Yaw measurement error'])
    
        legend( ...
            'Measured - Real', ...
            'Zero error', ...
            'Location','best')
    
    end
end

%% ================================================================
% 11C. REALE vs STIMATO DAL FILTRO DI KALMAN
%
% Confronto tra:
%   - stato reale
%   - stato stimato dal filtro di Kalman
%
% Solo alpha e beta.
%% ================================================================

figure( ...
    'Name','LQG comparison - Kalman estimate vs Real - Linear', ...
    'Color','w');

% ------------------------------------------------
% PITCH
% ------------------------------------------------
subplot(2,1,1)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tAlpha, ...
        rad2deg(D(i).alphaAbsV), ...
        'LineWidth',1.8, ...
        'DisplayName',[controllerNames{i} ' - Reale']);

    plot( ...
        D(i).tXhat, ...
        rad2deg(D(i).alphaHatAbsV), ...
        '--', ...
        'LineWidth',1.3, ...
        'DisplayName',[controllerNames{i} ' - Stimata']);

end

grid on
xlabel('Time [s]')
ylabel('\alpha [deg]')
title('Filtro di Kalman - Pitch: reale vs stimato')
legend('Location','best')


% ------------------------------------------------
% YAW
% ------------------------------------------------
subplot(2,1,2)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tBeta, ...
        rad2deg(D(i).betaAbsV), ...
        'LineWidth',1.8, ...
        'DisplayName',[controllerNames{i} ' - Reale']);

    plot( ...
        D(i).tXhat, ...
        rad2deg(D(i).betaHatAbsV), ...
        '--', ...
        'LineWidth',1.3, ...
        'DisplayName',[controllerNames{i} ' - Stimata']);

end

grid on
xlabel('Time [s]')
ylabel('\beta [deg]')
title('Filtro di Kalman - Yaw: reale vs stimato')
legend('Location','best')

%% ================================================================
% 11D. REALE vs STIMATO DAL FILTRO DI KALMAN - NON LINEARE
%% ================================================================

figure( ...
    'Name','LQG comparison - Kalman estimate vs Real - Nonlinear', ...
    'Color','w');

% ------------------------------------------------
% PITCH
% ------------------------------------------------
subplot(2,1,1)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tAlphaNL, ...
        rad2deg(D(i).alphaNLAbsV), ...
        'LineWidth',1.8, ...
        'DisplayName',[controllerNames{i} ' - Reale']);

    plot( ...
        D(i).tXhatNL, ...
        rad2deg(D(i).alphaHatNLAbsV), ...
        '--', ...
        'LineWidth',1.3, ...
        'DisplayName',[controllerNames{i} ' - Stimata']);

end

grid on
xlabel('Time [s]')
ylabel('\alpha [deg]')
title('Filtro di Kalman - Pitch: reale vs stimato - Non lineare')
legend('Location','best')


% ------------------------------------------------
% YAW
% ------------------------------------------------
subplot(2,1,2)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tBetaNL, ...
        rad2deg(D(i).betaNLAbsV), ...
        'LineWidth',1.8, ...
        'DisplayName',[controllerNames{i} ' - Reale']);

    plot( ...
        D(i).tXhatNL, ...
        rad2deg(D(i).betaHatNLAbsV), ...
        '--', ...
        'LineWidth',1.3, ...
        'DisplayName',[controllerNames{i} ' - Stimata']);

end

grid on
xlabel('Time [s]')
ylabel('\beta [deg]')
title('Filtro di Kalman - Yaw: reale vs stimato - Non lineare')
legend('Location','best')

%% ================================================================
% 11E. ERRORE DI STIMA DEL FILTRO DI KALMAN - LINEARE
%
% e_hat = x - x_hat
%% ================================================================

figure( ...
    'Name','LQG comparison - Kalman estimation error - Linear', ...
    'Color','w');

% ------------------------------------------------
% PITCH
% ------------------------------------------------
subplot(2,1,1)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tXhat, ...
        rad2deg(D(i).errEstAlpha), ...
        'LineWidth',1.3, ...
        'DisplayName',controllerNames{i});

end

yline(0,'k--','LineWidth',0.8)

grid on
xlabel('Time [s]')
ylabel('$\alpha-\hat{\alpha}$ [deg]', ...
       'Interpreter','latex')
title('Errore di stima Kalman - Pitch')
legend('Location','best')



% ------------------------------------------------
% YAW
% ------------------------------------------------
subplot(2,1,2)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tXhat, ...
        rad2deg(D(i).errEstBeta), ...
        'LineWidth',1.3, ...
        'DisplayName',controllerNames{i});

end

yline(0,'k--','LineWidth',0.8)

grid on
xlabel('Time [s]')
ylabel('$\beta-\hat{\beta}$ [deg]', ...
       'Interpreter','latex')
title('Errore di stima Kalman - Yaw')
legend('Location','best')


%% ================================================================
% 11F. ERRORE DI STIMA DEL FILTRO DI KALMAN - NON LINEARE
%% ================================================================

figure( ...
    'Name','LQG comparison - Kalman estimation error - Nonlinear', ...
    'Color','w');

% ------------------------------------------------
% PITCH
% ------------------------------------------------
subplot(2,1,1)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tXhatNL, ...
        rad2deg(D(i).errEstAlphaNL), ...
        'LineWidth',1.3, ...
        'DisplayName',controllerNames{i});

end

yline(0,'k--','LineWidth',0.8)

grid on
xlabel('Time [s]')
ylabel('$\alpha-\hat{\alpha}$ [deg]', ...
       'Interpreter','latex')
title('Errore di stima Kalman - Pitch - Non lineare')
legend('Location','best')


% ------------------------------------------------
% YAW
% ------------------------------------------------
subplot(2,1,2)
hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);

    plot( ...
        D(i).tXhatNL, ...
        rad2deg(D(i).errEstBetaNL), ...
        'LineWidth',1.3, ...
        'DisplayName',controllerNames{i});

end

yline(0,'k--','LineWidth',0.8)

grid on
xlabel('Time [s]')
ylabel('$\beta-\hat{\beta}$ [deg]', ...
       'Interpreter','latex')
title('Errore di stima Kalman - Yaw - Non lineare')
legend('Location','best')

%% ================================================================
% 11NL. MISURATO vs REALE - MODELLO NON LINEARE
%% ================================================================

if tipo_test == 2 || tipo_test == 3

    
    for j = 1:length(idx_confronto)
        i = idx_confronto(j);
    
        figure( ...
            'Name',[controllerNames{i} ' - Measured vs Real - Nonlinear'], ...
            'Color','w')
    
    
        % ------------------------------------------------------------
        % PITCH
        % ------------------------------------------------------------
    
        subplot(2,1,1)
    
        plot( ...
            D(i).tAlphaMeasNL, ...
            rad2deg(D(i).alphaMeasNLAbsV), ...
            '.', ...
            'MarkerSize',3)
        
        hold on
        
        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(D(i).alphaNLAbsV), ...
            'k-', ...
            'LineWidth',2)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('\alpha [deg]')
    
        title([controllerNames{i} ' - Pitch: real vs measured - Nonlinear'])
    
        legend( ...
            'Measured', ...
            'Real', ...
            'Location','best')
    
    
        % ------------------------------------------------------------
        % YAW
        % ------------------------------------------------------------
    
        subplot(2,1,2)
    
        plot( ...
            D(i).tBetaMeasNL, ...
            rad2deg(D(i).betaMeasNLAbsV), ...
            '.', ...
            'MarkerSize',3)
        
        hold on
        
        plot( ...
            D(i).tBetaNL, ...
            rad2deg(D(i).betaNLAbsV), ...
            'k-', ...
            'LineWidth',2)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('\beta [deg]')
    
        title([controllerNames{i} ' - Yaw: real vs measured - Nonlinear'])
    
        legend( ...
            'Measured', ...
            'Real', ...
            'Location','best')
    
    end

end

%% ================================================================
% 11BN. ERRORE DI MISURA - MODELLO NON LINEARE
%
% e_alpha = alpha_misurata - alpha_reale
% e_beta  = beta_misurata  - beta_reale
%% ================================================================

if tipo_test == 2 || tipo_test == 3

    
    for j = 1:length(idx_confronto)
        i = idx_confronto(j);
    
        % -------------------------------------------------------------
        % Errore di misura
        % -------------------------------------------------------------
    
        eAlphaMeasNL = ...
            D(i).alphaMeasNLV - D(i).alphaNLV;
    
        eBetaMeasNL = ...
            D(i).betaMeasNLV - D(i).betaNLV;
    
    
        figure( ...
            'Name',[controllerNames{i} ' - Measurement Error - Nonlinear'], ...
            'Color','w')
    
    
        % -------------------------------------------------------------
        % PITCH
        % -------------------------------------------------------------
    
        subplot(2,1,1)
    
        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(eAlphaMeasNL), ...
            'LineWidth',1.3)
    
        hold on
    
        yline(0,'k--','LineWidth',0.8)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('e_\alpha [deg]')
    
        title([controllerNames{i} ' - Pitch measurement error - Nonlinear'])
    
        legend( ...
            'Measured - Real', ...
            'Zero error', ...
            'Location','best')
    
    
        % -------------------------------------------------------------
        % YAW
        % -------------------------------------------------------------
    
        subplot(2,1,2)
    
        plot( ...
            D(i).tBetaNL, ...
            rad2deg(eBetaMeasNL), ...
            'LineWidth',1.3)
    
        hold on
    
        yline(0,'k--','LineWidth',0.8)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('e_\beta [deg]')
    
        title([controllerNames{i} ' - Yaw measurement error - Nonlinear'])
    
        legend( ...
            'Measured - Real', ...
            'Zero error', ...
            'Location','best')
    
    end


end


%% ================================================================
% 12. CONFRONTO TRACKING LINEARE - TUTTI I CONTROLLORI
%% ================================================================

figure( ...
    'Name','LQG comparison - Linear tracking', ...
    'Color','w')


subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tAlpha, ...
        rad2deg(D(i).alphaAbsV), ...
        'LineWidth',1.5)

end

yline( ...
    rad2deg(r_alpha_abs), ...
    'k:', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\alpha [deg]')

title('Confronto tracking lineare - Pitch')

legend( ...
    labelsConfronto{:}, ...
    'Reference', ...
    'Location','best')


subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tBeta, ...
        rad2deg(D(i).betaAbsV), ...
        'LineWidth',1.5)

end

yline( ...
    rad2deg(r_beta_abs), ...
    'k:', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\beta [deg]')

title('Confronto tracking lineare - Yaw')

legend( ...
    labelsConfronto{:}, ...
    'Reference', ...
    'Location','best')


%% ================================================================
% 13. CONFRONTO TRACKING NON LINEARE - TUTTI I CONTROLLORI
%% ================================================================

figure( ...
    'Name','LQG comparison - Nonlinear tracking', ...
    'Color','w')


subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tAlphaNL, ...
        rad2deg(D(i).alphaNLAbsV), ...
        'LineWidth',1.5)

end

yline( ...
    rad2deg(r_alpha_abs), ...
    'k:', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\alpha [deg]')

title('Confronto tracking non lineare - Pitch')

legend( ...
    labelsConfronto{:}, ...
    'Reference', ...
    'Location','best')


subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tBetaNL, ...
        rad2deg(D(i).betaNLAbsV), ...
        'LineWidth',1.5)

end

yline( ...
    rad2deg(r_beta_abs), ...
    'k:', ...
    'LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('\beta [deg]')

title('Confronto tracking non lineare - Yaw')

legend( ...
    labelsConfronto{:}, ...
    'Reference', ...
    'Location','best')


%% ================================================================
% 14. CONFRONTO ERRORE DI TRACKING - TUTTI I CONTROLLORI
%% ================================================================

figure( ...
    'Name','LQG comparison - Tracking error', ...
    'Color','w')


subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    errAlpha = r_alpha - D(i).alphaV;

    plot( ...
        D(i).tAlpha, ...
        rad2deg(errAlpha), ...
        'LineWidth',1.5)

end

yline(0,'k:','LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('e_\alpha [deg]')

title('Confronto errore di tracking - Modello lineare')

legend( ...
    labelsConfronto{:}, ...
    'Zero error', ...
    'Location','best')


subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    errBeta = r_beta - D(i).betaV;

    plot( ...
        D(i).tBeta, ...
        rad2deg(errBeta), ...
        'LineWidth',1.5)

end

yline(0,'k:','LineWidth',1)

grid on

xlabel('Time [s]')
ylabel('e_\beta [deg]')

title('Confronto errore di tracking Yaw - Modello lineare')

legend( ...
    labelsConfronto{:}, ...
    'Zero error', ...
    'Location','best')

%% ================================================================
% 14B. DISTURBANCE RECOVERY - TUTTI I CONTROLLORI
%      Confronto modello lineare e non lineare
%
%      Solo per tipo_test = 2 oppure 3
%% ================================================================

if tipo_test == 2 || tipo_test == 3

    t_dist = 10;

    % ============================================================
    % DISTURBANCE RECOVERY - MODELLO LINEARE
    % ============================================================

    figure( ...
        'Name','LQG comparison - Disturbance recovery - Linear', ...
        'Color','w')

    subplot(2,1,1)
    hold on

    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);

        errAlpha = r_alpha - D(i).alphaV;

        plot( ...
            D(i).tAlpha, ...
            rad2deg(errAlpha), ...
            'LineWidth',1.5)

    end

    xline(t_dist,'k--','LineWidth',1.2)
    yline(0,'k:','LineWidth',0.8)

    grid on
    xlim([8 20])

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Disturbance recovery - Pitch - Modello lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Disturbance at 10 s', ...
        'Zero error', ...
        'Location','best')


    subplot(2,1,2)
    hold on

    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);

        errBeta = r_beta - D(i).betaV;

        plot( ...
            D(i).tBeta, ...
            rad2deg(errBeta), ...
            'LineWidth',1.5)

    end

    xline(t_dist,'k--','LineWidth',1.2)
    yline(0,'k:','LineWidth',0.8)

    grid on
    xlim([8 20])

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Disturbance recovery - Yaw - Modello lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Disturbance at 10 s', ...
        'Zero error', ...
        'Location','best')


    % ============================================================
    % DISTURBANCE RECOVERY - MODELLO NON LINEARE
    % ============================================================

    figure( ...
        'Name','LQG comparison - Disturbance recovery - Nonlinear', ...
        'Color','w')

    subplot(2,1,1)
    hold on

    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);

        errAlphaNL = r_alpha - D(i).alphaNLV;

        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(errAlphaNL), ...
            'LineWidth',1.5)

    end

    xline(t_dist,'k--','LineWidth',1.2)
    yline(0,'k:','LineWidth',0.8)

    grid on
    xlim([8 20])

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Disturbance recovery - Pitch - Modello non lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Disturbance at 10 s', ...
        'Zero error', ...
        'Location','best')


    subplot(2,1,2)
    hold on

    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);

        errBetaNL = r_beta - D(i).betaNLV;

        plot( ...
            D(i).tBetaNL, ...
            rad2deg(errBetaNL), ...
            'LineWidth',1.5)

    end

    xline(t_dist,'k--','LineWidth',1.2)
    yline(0,'k:','LineWidth',0.8)

    grid on
    xlim([8 20])

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Disturbance recovery - Yaw - Modello non lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Disturbance at 10 s', ...
        'Zero error', ...
        'Location','best')

end

%% ================================================================
% 15. CONFRONTO LQG 2-DOF vs LQGI - TRACKING
%
% Confronto specifico tra:
%   - LQG 2-DOF
%   - LQGI
%
% Effettuato sia sul modello lineare sia sul modello non lineare.
% Solo per Tipo 1 e Tipo 2.
%% ================================================================

if tipo_test == 1 || tipo_test == 2

    idxLQG  = 2;
    idxLQGI = 3;


    % ============================================================
    % MODELLO LINEARE
    % ============================================================

    figure( ...
        'Name','LQG 2-DOF vs LQGI - Linear Tracking', ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(idxLQG).tAlpha, ...
        rad2deg(D(idxLQG).alphaAbsV), ...
        'LineWidth',1.5)
    
    hold on
    
    plot( ...
        D(idxLQGI).tAlpha, ...
        rad2deg(D(idxLQGI).alphaAbsV), ...
        'LineWidth',1.5)
    
    yline( ...
        rad2deg(r_alpha_abs), ...
        'k:', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\alpha [deg]')

    title('LQG 2-DOF vs LQGI - Pitch - Modello lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Reference', ...
        'Location','best')


    subplot(2,1,2)

    plot( ...
        D(idxLQG).tBeta, ...
        rad2deg(D(idxLQG).betaAbsV), ...
        'LineWidth',1.5)
    
    hold on
    
    plot( ...
        D(idxLQGI).tBeta, ...
        rad2deg(D(idxLQGI).betaAbsV), ...
        'LineWidth',1.5)
    yline( ...
        rad2deg(r_beta_abs), ...
        'k:', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\beta [deg]')

    title('LQG 2-DOF vs LQGI - Yaw - Modello lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Reference', ...
        'Location','best')


    % ============================================================
    % MODELLO NON LINEARE
    % ============================================================

    figure( ...
        'Name','LQG 2-DOF vs LQGI - Nonlinear Tracking', ...
        'Color','w')


    subplot(2,1,1)

    plot( ...
        D(idxLQG).tAlphaNL, ...
        rad2deg(D(idxLQG).alphaNLAbsV), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tAlphaNL, ...
        rad2deg(D(idxLQGI).alphaNLAbsV), ...
        'LineWidth',1.5)

    yline( ...
        rad2deg(r_alpha_abs), ...
        'k:', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\alpha [deg]')

    title('LQG 2-DOF vs LQGI - Pitch - Modello non lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Reference', ...
        'Location','best')


    subplot(2,1,2)

    plot( ...
        D(idxLQG).tBetaNL, ...
        rad2deg(D(idxLQG).betaNLAbsV), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tBetaNL, ...
        rad2deg(D(idxLQGI).betaNLAbsV), ...
        'LineWidth',1.5)

    yline( ...
        rad2deg(r_beta_abs), ...
        'k:', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\beta [deg]')

    title('LQG 2-DOF vs LQGI - Yaw - Modello non lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Reference', ...
        'Location','best')

end


%% ================================================================
% 16. CONFRONTO ERRORE LQG 2-DOF vs LQGI
%
% Sia sul modello lineare sia sul modello non lineare.
% Solo per Tipo 1 e Tipo 2.
%% ================================================================

if tipo_test == 1 || tipo_test == 2

    idxLQG  = 2;
    idxLQGI = 3;


    % ============================================================
    % MODELLO LINEARE
    % ============================================================

    figure( ...
        'Name','LQG 2-DOF vs LQGI - Linear Tracking Error', ...
        'Color','w')


    subplot(2,1,1)

    errAlphaLQG = ...
        r_alpha - D(idxLQG).alphaV;

    errAlphaLQGI = ...
        r_alpha - D(idxLQGI).alphaV;

    plot( ...
        D(idxLQG).tAlpha, ...
        rad2deg(errAlphaLQG), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tAlpha, ...
        rad2deg(errAlphaLQGI), ...
        'LineWidth',1.5)

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Errore di tracking - Pitch - Modello lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Zero error', ...
        'Location','best')


    subplot(2,1,2)

    errBetaLQG = ...
        r_beta - D(idxLQG).betaV;

    errBetaLQGI = ...
        r_beta - D(idxLQGI).betaV;

    plot( ...
        D(idxLQG).tBeta, ...
        rad2deg(errBetaLQG), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tBeta, ...
        rad2deg(errBetaLQGI), ...
        'LineWidth',1.5)

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Errore di tracking - Yaw - Modello lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Zero error', ...
        'Location','best')


    % ============================================================
    % MODELLO NON LINEARE
    % ============================================================

    figure( ...
        'Name','LQG 2-DOF vs LQGI - Nonlinear Tracking Error', ...
        'Color','w')


    subplot(2,1,1)

    errAlphaLQG_NL = ...
        r_alpha - D(idxLQG).alphaNLV;

    errAlphaLQGI_NL = ...
        r_alpha - D(idxLQGI).alphaNLV;

    plot( ...
        D(idxLQG).tAlphaNL, ...
        rad2deg(errAlphaLQG_NL), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tAlphaNL, ...
        rad2deg(errAlphaLQGI_NL), ...
        'LineWidth',1.5)

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Errore di tracking - Pitch - Modello non lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Zero error', ...
        'Location','best')


    subplot(2,1,2)

    errBetaLQG_NL = ...
        r_beta - D(idxLQG).betaNLV;

    errBetaLQGI_NL = ...
        r_beta - D(idxLQGI).betaNLV;

    plot( ...
        D(idxLQG).tBetaNL, ...
        rad2deg(errBetaLQG_NL), ...
        'LineWidth',1.5)

    hold on

    plot( ...
        D(idxLQGI).tBetaNL, ...
        rad2deg(errBetaLQGI_NL), ...
        'LineWidth',1.5)

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Errore di tracking - Yaw - Modello non lineare')

    legend( ...
        'LQG 2-DOF', ...
        'LQGI', ...
        'Zero error', ...
        'Location','best')

end

%% ================================================================
% 17. CONFRONTO SFORZO DI CONTROLLO
%     Modello lineare e non lineare
%% ================================================================

% ================================================================
% MODELLO LINEARE
% ================================================================

figure( ...
'Name','LQG comparison - Control effort - Linear', ...
'Color','w')

subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);
    
    
    plot( ...
        D(i).tUcmd, ...
        D(i).uCmdV(:,1), ...
        'LineWidth',1.4, ...
        'DisplayName',controllerNames{i})


end

grid on

xlabel('Time [s]')
ylabel('\Delta F_1 [N]')

title('Confronto comando di controllo - Main rotor - Modello lineare')

legend( ...
'Location','best')

subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);
    
    
    plot( ...
        D(i).tUcmd, ...
        D(i).uCmdV(:,2), ...
        'LineWidth',1.4, ...
        'DisplayName',controllerNames{i})


end

grid on

xlabel('Time [s]')
ylabel('\Delta F_2 [N]')

title('Confronto comando di controllo - Tail rotor - Modello lineare')

legend( ...
'Location','best')

% ================================================================
% MODELLO NON LINEARE
% ================================================================

figure( ...
'Name','LQG comparison - Control effort - Nonlinear', ...
'Color','w')

subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);
    
    
    plot( ...
        D(i).tUcmdNL, ...
        D(i).uCmdNLV(:,1), ...
        'LineWidth',1.4, ...
        'DisplayName',controllerNames{i})


end

grid on

xlabel('Time [s]')
ylabel('\Delta F_1 [N]')

title('Confronto comando di controllo - Main rotor - Modello non lineare')

legend( ...
'Location','best')

subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)
    i = idx_confronto(j);
    
    
    plot( ...
        D(i).tUcmdNL, ...
        D(i).uCmdNLV(:,2), ...
        'LineWidth',1.4, ...
        'DisplayName',controllerNames{i})


end

grid on

xlabel('Time [s]')
ylabel('\Delta F_2 [N]')

title('Confronto comando di controllo - Tail rotor - Modello non lineare')

legend( ...
'Location','best')


%% ================================================================
% 18. CONFRONTO USCITA REALE ATTUATORI
%% ================================================================

figure( ...
    'Name','LQG comparison - Actual actuator output', ...
    'Color','w')


subplot(2,1,1)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tDeltaF, ...
        D(i).deltaFV(:,1), ...
        'LineWidth',1.4)

end

grid on

xlabel('Time [s]')
ylabel('\Delta F_1 [N]')

title('Confronto uscita reale attuatore - Main rotor')

legend( ...
    labelsConfronto{:}, ...
    'Location','best')


subplot(2,1,2)

hold on

for j = 1:length(idx_confronto)      
    i = idx_confronto(j);

    plot( ...
        D(i).tDeltaF, ...
        D(i).deltaFV(:,2), ...
        'LineWidth',1.4)

end

grid on

xlabel('Time [s]')
ylabel('\Delta F_2 [N]')

title('Confronto uscita reale attuatore - Tail rotor')

legend( ...
    labelsConfronto{:}, ...
    'Location','best')


%% ================================================================
% 19. DISTURBI AERODINAMICI
%% ================================================================
if tipo_test == 2 || tipo_test == 3
    if hasAero
    
        figure( ...
            'Name','Aerodynamic disturbances', ...
            'Color','w')
    
    
        subplot(2,1,1)
    
        plot( ...
            t_aero, ...
            dist_aero_alpha, ...
            'LineWidth',1.5)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('d_\alpha [N m]')
    
        title('Pitch aerodynamic disturbance')
    
    
        subplot(2,1,2)
    
        plot( ...
            t_aero, ...
            dist_aero_beta, ...
            'LineWidth',1.5)
    
        grid on
    
        xlabel('Time [s]')
        ylabel('d_\beta [N m]')
    
        title('Yaw aerodynamic disturbance')
    
    end
end

if tipo_test == 2 || tipo_test == 3

    %% ================================================================
    % 20. CONFRONTO MISURATO vs REALE - TUTTI I CONTROLLORI
    %% ================================================================
    
    figure( ...
        'Name','LQG comparison - Measured vs Real', ...
        'Color','w');
    
    % ------------------------------------------------
    % PITCH
    % ------------------------------------------------
    subplot(2,1,1)
    hold on
    
    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);
    
        plot( ...
            D(i).tAlpha, ...
            rad2deg(D(i).alphaAbsV), ...
            '-', ...
            'LineWidth',1.8, ...
            'DisplayName',[controllerNames{i} ' - Reale']);
    
        plot( ...
            D(i).tAlphaMeas, ...
            rad2deg(D(i).alphaMeasAbsV), ...
            '.', ...
            'MarkerSize',4, ...
            'DisplayName',[controllerNames{i} ' - Misurata']);
    
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\alpha [deg]')
    title('Pitch: misurata vs reale')
    legend('Location','best')
    
    % ------------------------------------------------
    % YAW
    % ------------------------------------------------
    subplot(2,1,2)
    hold on
    
    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);
    
        plot( ...
            D(i).tBeta, ...
            rad2deg(D(i).betaAbsV), ...
            '-', ...
            'LineWidth',1.8, ...
            'DisplayName',[controllerNames{i} ' - Reale']);
    
        plot( ...
            D(i).tBetaMeas, ...
            rad2deg(D(i).betaMeasAbsV), ...
            '.', ...
            'MarkerSize',4, ...
            'DisplayName',[controllerNames{i} ' - Misurata']);
    
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\beta [deg]')
    title('Yaw: misurata vs reale')
    legend('Location','best')
    
    %% ================================================================
    % CONFRONTO ERRORE DI MISURA - TUTTI I CONTROLLORI
    %% ================================================================
    
    figure( ...
        'Name','LQG measurement error comparison', ...
        'Color','w');
    
    subplot(2,1,1)
    hold on
    
    for j = 1:length(idx_confronto)     
        i = idx_confronto(j);
        e_alpha = D(i).alphaMeasV - D(i).alphaV;
    
        plot( ...
            D(i).tAlpha, ...
            rad2deg(e_alpha), ...
            'LineWidth',1.2, ...
            'DisplayName',controllerNames{i});
    end
    
    yline(0,'k--')
    grid on
    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')
    title('Errore di misura - Pitch')
    legend('Location','best')
    
    subplot(2,1,2)
    hold on
    
    for j = 1:length(idx_confronto)      
        i = idx_confronto(j);
        e_beta = D(i).betaMeasV - D(i).betaV;
    
        plot( ...
            D(i).tBeta, ...
            rad2deg(e_beta), ...
            'LineWidth',1.2, ...
            'DisplayName',controllerNames{i});
    end
    
    yline(0,'k--')
    grid on
    xlabel('Time [s]')
    ylabel('e_\beta [deg]')
    title('Errore di misura - Yaw')
    legend('Location','best')

end

%% ================================================================
% 20NL. CONFRONTI GLOBALI - MODELLO NON LINEARE
%
%   20NL-A  -> Measured vs Real
%   20NL-B  -> Tracking Error
%   20NL-C  -> Measurement Error
%
% Solo per Tipo 2 e Tipo 3
%% ================================================================

if tipo_test == 2 || tipo_test == 3

    %% ============================================================
    % 20NL-A. CONFRONTO MISURATO vs REALE - NON LINEARE
    %% ============================================================

    figure( ...
        'Name','LQG comparison - Measured vs Real - Nonlinear', ...
        'Color','w');

    % ------------------------------------------------------------
    % PITCH
    % ------------------------------------------------------------

    subplot(2,1,1)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(D(i).alphaNLAbsV), ...
            '-', ...
            'LineWidth',1.8, ...
            'DisplayName',[controllerNames{i} ' - Reale']);

        plot( ...
            D(i).tAlphaMeasNL, ...
            rad2deg(D(i).alphaMeasNLAbsV), ...
            '.', ...
            'MarkerSize',4, ...
            'DisplayName',[controllerNames{i} ' - Misurata']);

    end

    grid on

    xlabel('Time [s]')
    ylabel('\alpha [deg]')

    title('Pitch: misurata vs reale - Modello non lineare')

    legend('Location','best')


    % ------------------------------------------------------------
    % YAW
    % ------------------------------------------------------------

    subplot(2,1,2)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        plot( ...
            D(i).tBetaNL, ...
            rad2deg(D(i).betaNLAbsV), ...
            '-', ...
            'LineWidth',1.8, ...
            'DisplayName',[controllerNames{i} ' - Reale']);

        plot( ...
            D(i).tBetaMeasNL, ...
            rad2deg(D(i).betaMeasNLAbsV), ...
            '.', ...
            'MarkerSize',4, ...
            'DisplayName',[controllerNames{i} ' - Misurata']);

    end

    grid on

    xlabel('Time [s]')
    ylabel('\beta [deg]')

    title('Yaw: misurata vs reale - Modello non lineare')

    legend('Location','best')


    %% ============================================================
    % 20NL-B. CONFRONTO ERRORE DI TRACKING - NON LINEARE
    %% ============================================================

    figure( ...
        'Name','LQG comparison - Tracking error - Nonlinear', ...
        'Color','w');

    % ------------------------------------------------------------
    % PITCH
    % ------------------------------------------------------------

    subplot(2,1,1)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        errAlphaNL = ...
            r_alpha - D(i).alphaNLV;

        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(errAlphaNL), ...
            'LineWidth',1.5, ...
            'DisplayName',controllerNames{i});

    end

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Confronto errore di tracking - Pitch - Modello non lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Zero error', ...
        'Location','best')


    % ------------------------------------------------------------
    % YAW
    % ------------------------------------------------------------

    subplot(2,1,2)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        errBetaNL = ...
            r_beta - D(i).betaNLV;

        plot( ...
            D(i).tBetaNL, ...
            rad2deg(errBetaNL), ...
            'LineWidth',1.5, ...
            'DisplayName',controllerNames{i});

    end

    yline(0,'k:','LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Confronto errore di tracking - Yaw - Modello non lineare')

    legend( ...
        labelsConfronto{:}, ...
        'Zero error', ...
        'Location','best')


    %% ============================================================
    % 20NL-C. CONFRONTO ERRORE DI MISURA - NON LINEARE
    %% ============================================================

    figure( ...
        'Name','LQG measurement error comparison - Nonlinear', ...
        'Color','w');

    % ------------------------------------------------------------
    % PITCH
    % ------------------------------------------------------------

    subplot(2,1,1)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        e_alphaNL = ...
            D(i).alphaMeasNLV - D(i).alphaNLV;

        plot( ...
            D(i).tAlphaNL, ...
            rad2deg(e_alphaNL), ...
            'LineWidth',1.2, ...
            'DisplayName',controllerNames{i});

    end

    yline(0,'k--','LineWidth',0.8)

    grid on

    xlabel('Time [s]')
    ylabel('e_\alpha [deg]')

    title('Errore di misura - Pitch - Modello non lineare')

    legend('Location','best')


    % ------------------------------------------------------------
    % YAW
    % ------------------------------------------------------------

    subplot(2,1,2)
    hold on

    for j = 1:length(idx_confronto)

        i = idx_confronto(j);

        e_betaNL = ...
            D(i).betaMeasNLV - D(i).betaNLV;

        plot( ...
            D(i).tBetaNL, ...
            rad2deg(e_betaNL), ...
            'LineWidth',1.2, ...
            'DisplayName',controllerNames{i});

    end

    yline(0,'k--','LineWidth',0.8)

    grid on

    xlabel('Time [s]')
    ylabel('e_\beta [deg]')

    title('Errore di misura - Yaw - Modello non lineare')

    legend('Location','best')

end

%% ================================================================
% 21. TABELLA RIASSUNTIVA
%% ================================================================

NomeControllore = controllerNames(idx_confronto).';

Nconfronto = length(idx_confronto);

RMS_alpha_lin = zeros(Nconfronto,1);
RMS_beta_lin  = zeros(Nconfronto,1);

RMS_alpha_nl = zeros(Nconfronto,1);
RMS_beta_nl  = zeros(Nconfronto,1);

Peak_alpha_lin = zeros(Nconfronto,1);
Peak_beta_lin  = zeros(Nconfronto,1);

Peak_alpha_nl = zeros(Nconfronto,1);
Peak_beta_nl  = zeros(Nconfronto,1);

Final_alpha_lin = zeros(Nconfronto,1);
Final_beta_lin  = zeros(Nconfronto,1);

Final_alpha_nl = zeros(Nconfronto,1);
Final_beta_nl  = zeros(Nconfronto,1);

% Metriche di transitorio - modello lineare
Overshoot_alpha_lin = zeros(Nconfronto,1);
Rise_alpha_lin      = zeros(Nconfronto,1);
Settling_alpha_lin  = zeros(Nconfronto,1);

% Metriche di transitorio - modello non lineare
Overshoot_alpha_nl = zeros(Nconfronto,1);
Rise_alpha_nl      = zeros(Nconfronto,1);
Settling_alpha_nl  = zeros(Nconfronto,1);


for j = 1:length(idx_confronto)     
    i = idx_confronto(j);

    RMS_alpha_lin(j) = D(i).ris_lin.RMS_alpha;
    RMS_beta_lin(j)  = D(i).ris_lin.RMS_beta;

    RMS_alpha_nl(j) = D(i).ris_nl.RMS_alpha;
    RMS_beta_nl(j)  = D(i).ris_nl.RMS_beta;

    Peak_alpha_lin(j) = D(i).ris_lin.peak_alpha;
    Peak_beta_lin(j)  = D(i).ris_lin.peak_beta;

    Peak_alpha_nl(j) = D(i).ris_nl.peak_alpha;
    Peak_beta_nl(j)  = D(i).ris_nl.peak_beta;

    Final_alpha_lin(j) = D(i).ris_lin.final_error_alpha;
    Final_beta_lin(j)  = D(i).ris_lin.final_error_beta;

    Final_alpha_nl(j) = D(i).ris_nl.final_error_alpha;
    Final_beta_nl(j) = D(i).ris_nl.final_error_beta;

    % Metriche transitorie lineari
    Overshoot_alpha_lin(j) = D(i).ris_lin.overshoot_alpha;
    Rise_alpha_lin(j)      = D(i).ris_lin.rise_time_alpha;
    Settling_alpha_lin(j)  = D(i).ris_lin.settling_time_alpha;
    
    % Metriche transitorie non lineari
    Overshoot_alpha_nl(j) = D(i).ris_nl.overshoot_alpha;
    Rise_alpha_nl(j)      = D(i).ris_nl.rise_time_alpha;
    Settling_alpha_nl(j)  = D(i).ris_nl.settling_time_alpha;

end


Tabella_LQG = table( ...
    NomeControllore, ...
    RMS_alpha_lin, ...
    RMS_beta_lin, ...
    RMS_alpha_nl, ...
    RMS_beta_nl, ...
    Peak_alpha_lin, ...
    Peak_beta_lin, ...
    Peak_alpha_nl, ...
    Peak_beta_nl, ...
    Final_alpha_lin, ...
    Final_beta_lin, ...
    Final_alpha_nl, ...
    Final_beta_nl, ...
    Rise_alpha_lin, ...
    Overshoot_alpha_lin, ...
    Settling_alpha_lin, ...
    Rise_alpha_nl, ...
    Overshoot_alpha_nl, ...
    Settling_alpha_nl);

fprintf('\n\n')
fprintf('############################################################\n')
fprintf('# TABELLA RIASSUNTIVA LQG / LQGI\n')
fprintf('############################################################\n')

disp(Tabella_LQG)


%% ================================================================
% 22. RMS ERRORE DI LINEARIZZAZIONE
%% ================================================================

RMS_lin_alpha = zeros(Nconfronto,1);
RMS_lin_beta  = zeros(Nconfronto,1);

Peak_lin_alpha = zeros(Nconfronto,1);
Peak_lin_beta  = zeros(Nconfronto,1);


for j = 1:length(idx_confronto)     
    i = idx_confronto(j);

    RMS_lin_alpha(j) = sqrt( ...
        mean(D(i).errAlpha.^2));

    RMS_lin_beta(j) = sqrt( ...
        mean(D(i).errBeta.^2));

    Peak_lin_alpha(j) = max(abs(D(i).errAlpha));

    Peak_lin_beta(j) = max(abs(D(i).errBeta));

end


Tabella_linearizzazione = table( ...
    NomeControllore, ...
    RMS_lin_alpha, ...
    RMS_lin_beta, ...
    Peak_lin_alpha, ...
    Peak_lin_beta);


fprintf('\n')
fprintf('============================================================\n')
fprintf(' ERRORE DI LINEARIZZAZIONE\n')
fprintf('============================================================\n')

disp(Tabella_linearizzazione)

%% ================================================================
% 23. METRICHE ERRORE DI STIMA KALMAN
%
% e_hat = x - x_hat
%% ================================================================

RMS_est_alpha_lin = zeros(Nconfronto,1);
RMS_est_beta_lin  = zeros(Nconfronto,1);

Peak_est_alpha_lin = zeros(Nconfronto,1);
Peak_est_beta_lin  = zeros(Nconfronto,1);

RMS_est_alpha_nl = zeros(Nconfronto,1);
RMS_est_beta_nl  = zeros(Nconfronto,1);

Peak_est_alpha_nl = zeros(Nconfronto,1);
Peak_est_beta_nl  = zeros(Nconfronto,1);


for j = 1:length(idx_confronto)

    i = idx_confronto(j);

    % -----------------------------------------
    % Modello lineare
    % -----------------------------------------

    RMS_est_alpha_lin(j) = rms(D(i).errEstAlpha);
    RMS_est_beta_lin(j)  = rms(D(i).errEstBeta);

    Peak_est_alpha_lin(j) = max(abs(D(i).errEstAlpha));
    Peak_est_beta_lin(j)  = max(abs(D(i).errEstBeta));


    % -----------------------------------------
    % Modello non lineare
    % -----------------------------------------

    RMS_est_alpha_nl(j) = rms(D(i).errEstAlphaNL);
    RMS_est_beta_nl(j)  = rms(D(i).errEstBetaNL);

    Peak_est_alpha_nl(j) = max(abs(D(i).errEstAlphaNL));
    Peak_est_beta_nl(j)  = max(abs(D(i).errEstBetaNL));

end

% ================================================================
% DIAGNOSTICA ERRORE KALMAN
% ================================================================

for j = 1:length(idx_confronto)

    i = idx_confronto(j);

    fprintf('\n===============================================\n')
    fprintf('%s - DIAGNOSTICA KALMAN\n',controllerNames{i})
    fprintf('===============================================\n')

    % LINEARE
    meanAlpha = mean(D(i).errEstAlpha);
    meanBeta  = mean(D(i).errEstBeta);

    stdAlpha = std(D(i).errEstAlpha);
    stdBeta  = std(D(i).errEstBeta);

    fprintf('\nMODELLO LINEARE\n')
    fprintf('Mean error alpha = %.6f rad (%.4f deg)\n', ...
        meanAlpha,rad2deg(meanAlpha))
    fprintf('Mean error beta  = %.6f rad (%.4f deg)\n', ...
        meanBeta,rad2deg(meanBeta))

    fprintf('STD error alpha  = %.6f rad (%.4f deg)\n', ...
        stdAlpha,rad2deg(stdAlpha))
    fprintf('STD error beta   = %.6f rad (%.4f deg)\n', ...
        stdBeta,rad2deg(stdBeta))


    % NON LINEARE
    meanAlphaNL = mean(D(i).errEstAlphaNL);
    meanBetaNL  = mean(D(i).errEstBetaNL);

    stdAlphaNL = std(D(i).errEstAlphaNL);
    stdBetaNL  = std(D(i).errEstBetaNL);

    fprintf('\nMODELLO NON LINEARE\n')
    fprintf('Mean error alpha = %.6f rad (%.4f deg)\n', ...
        meanAlphaNL,rad2deg(meanAlphaNL))
    fprintf('Mean error beta  = %.6f rad (%.4f deg)\n', ...
        meanBetaNL,rad2deg(meanBetaNL))

    fprintf('STD error alpha  = %.6f rad (%.4f deg)\n', ...
        stdAlphaNL,rad2deg(stdAlphaNL))
    fprintf('STD error beta   = %.6f rad (%.4f deg)\n', ...
        stdBetaNL,rad2deg(stdBetaNL))

end

%% ================================================================
% ANALISI ERRORE KALMAN PRIMA E DOPO IL DISTURBO
% ================================================================

t_dist = 10;

for j = 1:length(idx_confronto)

    i = idx_confronto(j);

    t = D(i).tXhat;

    % ------------------------------------------------------------
    % PRE-DISTURBO
    % ------------------------------------------------------------

    idx_pre = t < t_dist;

    rms_alpha_pre = rms(D(i).errEstAlpha(idx_pre));
    rms_beta_pre  = rms(D(i).errEstBeta(idx_pre));

    peak_alpha_pre = max(abs(D(i).errEstAlpha(idx_pre)));
    peak_beta_pre  = max(abs(D(i).errEstBeta(idx_pre)));

    mean_alpha_pre = mean(D(i).errEstAlpha(idx_pre));
    mean_beta_pre  = mean(D(i).errEstBeta(idx_pre));


    % ------------------------------------------------------------
    % DOPO IL DISTURBO
    % ------------------------------------------------------------

    idx_post = t >= t_dist;

    rms_alpha_post = rms(D(i).errEstAlpha(idx_post));
    rms_beta_post  = rms(D(i).errEstBeta(idx_post));

    peak_alpha_post = max(abs(D(i).errEstAlpha(idx_post)));
    peak_beta_post  = max(abs(D(i).errEstBeta(idx_post)));

    mean_alpha_post = mean(D(i).errEstAlpha(idx_post));
    mean_beta_post  = mean(D(i).errEstBeta(idx_post));


    fprintf('\n')
    fprintf('============================================================\n')
    fprintf('%s - ERRORE KALMAN PRE/POST DISTURBO\n',controllerNames{i})
    fprintf('============================================================\n')

    fprintf('\nPRE-DISTURBO (t < %.1f s)\n',t_dist)

    fprintf('RMS alpha = %.6f rad (%.3f deg)\n', ...
        rms_alpha_pre,rad2deg(rms_alpha_pre));

    fprintf('RMS beta  = %.6f rad (%.3f deg)\n', ...
        rms_beta_pre,rad2deg(rms_beta_pre));

    fprintf('Peak alpha = %.6f rad (%.3f deg)\n', ...
        peak_alpha_pre,rad2deg(peak_alpha_pre));

    fprintf('Peak beta  = %.6f rad (%.3f deg)\n', ...
        peak_beta_pre,rad2deg(peak_beta_pre));

    fprintf('Mean alpha = %.6f rad (%.3f deg)\n', ...
        mean_alpha_pre,rad2deg(mean_alpha_pre));

    fprintf('Mean beta  = %.6f rad (%.3f deg)\n', ...
        mean_beta_pre,rad2deg(mean_beta_pre));


    fprintf('\nPOST-DISTURBO (t >= %.1f s)\n',t_dist)

    fprintf('RMS alpha = %.6f rad (%.3f deg)\n', ...
        rms_alpha_post,rad2deg(rms_alpha_post));

    fprintf('RMS beta  = %.6f rad (%.3f deg)\n', ...
        rms_beta_post,rad2deg(rms_beta_post));

    fprintf('Peak alpha = %.6f rad (%.3f deg)\n', ...
        peak_alpha_post,rad2deg(peak_alpha_post));

    fprintf('Peak beta  = %.6f rad (%.3f deg)\n', ...
        peak_beta_post,rad2deg(peak_beta_post));

    fprintf('Mean alpha = %.6f rad (%.3f deg)\n', ...
        mean_alpha_post,rad2deg(mean_alpha_post));

    fprintf('Mean beta  = %.6f rad (%.3f deg)\n', ...
        mean_beta_post,rad2deg(mean_beta_post));

end


Tabella_Kalman = table( ...
    NomeControllore, ...
    RMS_est_alpha_lin, ...
    RMS_est_beta_lin, ...
    Peak_est_alpha_lin, ...
    Peak_est_beta_lin, ...
    RMS_est_alpha_nl, ...
    RMS_est_beta_nl, ...
    Peak_est_alpha_nl, ...
    Peak_est_beta_nl);


fprintf('\n\n')
fprintf('############################################################\n')
fprintf('# METRICHE ERRORE DI STIMA - FILTRO DI KALMAN\n')
fprintf('############################################################\n')

disp(Tabella_Kalman)


%% ================================================================
% SALVATAGGIO AUTOMATICO DI TUTTE LE FIGURE
%
% Struttura:
%
% figure_LQG/
%
%   Tipo1/
%       LQG1/
%       LQG2/
%       LQI/
%       Confronti/
%
%   Tipo2/
%       LQG1/
%       LQG2/
%       LQI/
%       Confronti/
%
%   Tipo3/
%       LQG1/
%       LQG2/
%       LQI/
%       Confronti/
%
%% ================================================================

fprintf('\n')
fprintf('============================================================\n')
fprintf(' SALVATAGGIO FIGURE - TIPO %d\n',tipo_test)
fprintf('============================================================\n')


%% ================================================================
% 1. CARTELLE
%% ================================================================

cartellaBase = fullfile(pwd,'figure_LQG');

cartellaTipo = fullfile( ...
    cartellaBase, ...
    sprintf('Tipo%d',tipo_test));

cartellaLQG1 = fullfile(cartellaTipo,'LQG1');
cartellaLQG2 = fullfile(cartellaTipo,'LQG2');
cartellaLQI  = fullfile(cartellaTipo,'LQI');

cartellaConfronti = fullfile( ...
    cartellaTipo,'Confronti');


cartelle = { ...
    cartellaLQG1, ...
    cartellaLQG2, ...
    cartellaLQI, ...
    cartellaConfronti};


for k = 1:numel(cartelle)

    if ~exist(cartelle{k},'dir')
        mkdir(cartelle{k});
    end

end


%% ================================================================
% 2. TUTTE LE FIGURE APERTE
%% ================================================================

figures = findall(0,'Type','figure');

fprintf('Numero di figure trovate: %d\n\n',length(figures));


%% ================================================================
% 3. SALVATAGGIO
%% ================================================================

for k = 1:length(figures)

    fig = figures(k);

    nome = get(fig,'Name');

    if isempty(nome)
        nome = sprintf('Figura_%02d',k);
    end

    fprintf('\n[%d/%d] Sto salvando: %s\n',k,length(figures),nome);
    drawnow;


    % ------------------------------------------------------------
    % Nome file valido
    % ------------------------------------------------------------

    nomeFile = regexprep( ...
        nome, ...
        '[<>:"/\\|?*]', ...
        '_');

    nomeFile = strrep(nomeFile,' ','_');
    nomeFile = strrep(nomeFile,'-','_');

    nomeLower = lower(nome);


    % ------------------------------------------------------------
    % Destinazione
    % ------------------------------------------------------------

    isConfronto = ...
        contains(nomeLower,'comparison') || ...
        contains(nomeLower,'confronto') || ...
        contains(nomeLower,'vs lqg') || ...
        contains(nomeLower,'vs lqi') || ...
        contains(nomeLower,'vs lqgi') || ...
        contains(nomeLower,'vs_lqg') || ...
        contains(nomeLower,'vs_lqi') || ...
        contains(nomeLower,'vs_lqgi');

    if isConfronto

        cartellaDest = cartellaConfronti;

    elseif contains(nomeLower,'lqg 1-dof') || ...
           contains(nomeLower,'lqg1') || ...
           contains(nomeLower,'1dof')

        cartellaDest = cartellaLQG1;

    elseif contains(nomeLower,'lqg 2-dof') || ...
           contains(nomeLower,'lqg2') || ...
           contains(nomeLower,'2dof')

        cartellaDest = cartellaLQG2;

    elseif contains(nomeLower,'lqgi') || ...
           contains(nomeLower,'lqi')

        cartellaDest = cartellaLQI;

    else

        cartellaDest = cartellaConfronti;

    end


    % ------------------------------------------------------------
    % PNG
    % ------------------------------------------------------------

    fprintf('  -> exportgraphics...\n');
    drawnow;

    try
        exportgraphics( ...
            fig, ...
            fullfile(cartellaDest,[nomeFile '.png']), ...
            'Resolution',300);

        fprintf('  -> PNG salvato\n');
        drawnow;

    catch ME

        warning( ...
            'Impossibile salvare PNG per "%s": %s', ...
            nome,ME.message);

    end


    % ------------------------------------------------------------
    % FIG
    % ------------------------------------------------------------

    fprintf('  -> savefig...\n');
    drawnow;

    try

        savefig( ...
            fig, ...
            fullfile(cartellaDest,[nomeFile '.fig']));

        fprintf('  -> FIG salvato\n');
        drawnow;

    catch ME

        warning( ...
            'Impossibile salvare FIG per "%s": %s', ...
            nome,ME.message);

    end

    fprintf('Salvata: %s\n', ...
        fullfile(cartellaDest,[nomeFile '.png']));

end


fprintf('\n')
fprintf('Figure salvate in:\n%s\n',cartellaTipo)

fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('############################################################\n')
fprintf('# ANALISI LQG / LQGI COMPLETATA\n')
fprintf('############################################################\n')