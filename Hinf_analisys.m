clc

%% CONFIGURAZIONE

angleUnit = 'rad';       % 'rad' oppure 'deg'
tStepRef = 1;            % Istante del gradino alpha [s]
alphaRefDeg = 3;         % Riferimento alpha [deg]

timeVariable = '';       % Esempio: 'tout', se necessario

makeVisibleFigures = true;
outputDir = fullfile('images','Hinf','Confronti');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

controllers = {'mixsyn','hinfsyn','hinfstruct'};
labels = {'mixsyn','hinfsyn','hinfstruct'};
colors = lines(numel(controllers));

%% 1. TRACKING - MODELLO LINEARE

fig = figure('Name','Hinf - Linear tracking','Color','w', ...
    'Visible',onOff(makeVisibleFigures));
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plotTracking(controllers,labels,colors,'', ...
    angleUnit,tStepRef,alphaRefDeg,timeVariable);

saveFigure(fig,outputDir,'Hinf_comparison___Linear_tracking.png');

%% 2. TRACKING - MODELLO NON LINEARE

fig = figure('Name','Hinf - Nonlinear tracking','Color','w', ...
    'Visible',onOff(makeVisibleFigures));
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plotTracking(controllers,labels,colors,'_nl', ...
    angleUnit,tStepRef,alphaRefDeg,timeVariable);

saveFigure(fig,outputDir,'Hinf_comparison___Nonlinear_tracking.png');

%% 3. SFORZO DI CONTROLLO - MODELLO NON LINEARE

fig = figure('Name','Hinf - Control effort nonlinear','Color','w', ...
    'Visible',onOff(makeVisibleFigures));
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plotCommands(controllers,labels,colors,'_nl',timeVariable);

saveFigure(fig,outputDir, ...
    'Hinf_comparison___Control_effort___Nonlinear.png');

%% 4. ERRORE DI MISURA - MODELLO NON LINEARE

fig = figure('Name','Hinf - Measurement error nonlinear','Color','w', ...
    'Visible',onOff(makeVisibleFigures));
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plotMeasurementError(controllers,labels,colors, ...
    angleUnit,timeVariable);

saveFigure(fig,outputDir, ...
    'Hinf_measurement_error_comparison___Nonlinear.png');

%% 5. ERRORE DI LINEARIZZAZIONE

fig = figure('Name','Hinf - Linearization error','Color','w', ...
    'Visible',onOff(makeVisibleFigures));
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plotLinearizationError(controllers,labels,colors, ...
    angleUnit,timeVariable);

saveFigure(fig,outputDir,'Hinf_comparison___Linearization_Error.png');

%% 6. METRICHE DI TRACKING

metrics = buildTrackingMetrics(controllers,angleUnit, ...
    tStepRef,alphaRefDeg,timeVariable);

csvPath = fullfile(outputDir,'Hinf_tracking_metrics.csv');
writetable(metrics,csvPath);

fprintf('\nGenerazione completata.\nCartella: %s\n',outputDir);
disp(metrics);



%% ========================================================================
% GRAFICI E VERIFICA NOMINALE H-INFINITY
% ========================================================================

fprintf('\n=== VERIFICA NOMINALE H-INFINITY ===\n');

freqDir = fullfile(pwd,'images','Hinf');

if ~exist(freqDir,'dir')
    mkdir(freqDir);
end

% Controllori NORMALIZZATI: devono essere coerenti con G_scaled
controllers = {
    K_mix_scaled
    K_hinfsyn_scaled
    K_pidcomp_scaled
};

names = {
    'mixsyn'
    'hinfsyn'
    'hinfstruct'
};

w = logspace(-2,3,1000);
I = eye(size(G_scaled,1));

% Tabella risultati
results = zeros(3,4);

for k = 1:numel(controllers)

    K = controllers{k};

    % ------------------------------------------------------------
    % SISTEMI IN RETROAZIONE
    % ------------------------------------------------------------

    L = G_scaled*K;

    S = feedback(I,L);
    T = I-S;
    KS = K*S;

    % ------------------------------------------------------------
    % FUNZIONI PESATE
    % ------------------------------------------------------------

    WS_S = minreal(WS*S);
    WU_KS = minreal(WU*KS);
    WT_T = minreal(WT*T);

    % ------------------------------------------------------------
    % VERIFICA DELLE NORME H-INFINITY
    % ------------------------------------------------------------

    gamma_WS = hinfnorm(WS_S);
    gamma_WU = hinfnorm(WU_KS);
    gamma_WT = hinfnorm(WT_T);

    CL_weighted = [
        WS_S
        WU_KS
        WT_T
    ];

    gamma_total = hinfnorm(minreal(CL_weighted));

    results(k,:) = [
        gamma_WS
        gamma_WU
        gamma_WT
        gamma_total
    ];

    % ------------------------------------------------------------
    % GRAFICI DEI VALORI SINGOLARI PESATI
    % ------------------------------------------------------------

    weightedSystems = {
        WS_S
        WU_KS
        WT_T
    };

    plotTitles = {
        'Sensibilita'' pesata W_S S'
        'Sforzo di controllo pesato W_U KS'
        'Sensibilita'' complementare pesata W_T T'
    };

    fig = figure('Color','w', ...
        'Name',['Verifica Hinf - ' names{k}]);

    tiledlayout(3,1,'TileSpacing','compact');

    for j = 1:3

        nexttile;

        sig = maxSigmaDb(weightedSystems{j},w);

        semilogx(w,sig,'b','LineWidth',1.5);
        hold on;

        yline(0,'r--','LineWidth',1.2);

        grid on;
        xlim([w(1) w(end)]);
        ylim([-60 40]);

        ylabel('Modulo [dB]');
        title(plotTitles{j});

        legend('Valore singolare massimo', ...
               'Limite unitario','Location','best');

    end

    xlabel('Frequenza [rad/s]');

    exportgraphics(fig, ...
        fullfile(freqDir,['Sigma_' names{k} '_weighted.png']), ...
        'Resolution',300);

end

% ------------------------------------------------------------
% REPORT FINALE
% ------------------------------------------------------------

T_results = array2table(results, ...
    'VariableNames', ...
    {'Norm_WS_S','Norm_WU_KS','Norm_WT_T','Gamma_tot'}, ...
    'RowNames',names);

disp(' ');
disp('=== VERIFICA NOMINALE H-INFINITY ===');
disp(T_results);

fprintf('\nVerifica dei vincoli:\n');

for k = 1:3

    fprintf('\n%s:\n',names{k});

    if results(k,4) < 1
        fprintf('  TUTTI I VINCOLI RISPETTATI\n');
    else
        fprintf('  VINCOLO H-INFINITY NON RISPETTATO\n');
    end

end

writetable(T_results, ...
    fullfile(freqDir,'Hinf_frequency_verification.csv'), ...
    'WriteRowNames',true);

fprintf('\nGrafici e risultati salvati in: %s\n',freqDir);



%% ========================================================================
% ERRORE DI LINEARIZZAZIONE: MODELLO LINEARE vs NON LINEARE
% ========================================================================

controllers = {'mixsyn / hinfsyn'; 'hinfstruct'};

% Ogni riga: risposta lineare, risposta non lineare
pairs = {
    alpha_mixsyn,       alpha_mixsyn_nl, ...
    beta_mixsyn,        beta_mixsyn_nl;
    alpha_hinfstruct,   alpha_hinfstruct_nl, ...
    beta_hinfstruct,    beta_hinfstruct_nl
    };

% Colonne: RMS alpha, RMS beta, Peak alpha, Peak beta
metrics = zeros(2,4);

for i = 1:2

    % Alpha
    [t_lin, y_lin, ok_lin] = getSignal(pairs{i,1});
    [t_nl,  y_nl,  ok_nl]  = getSignal(pairs{i,2});

    [e_alpha, valid] = linearizationError(t_lin, y_lin, t_nl, y_nl);

    metrics(i,1) = sqrt(mean(e_alpha(valid).^2));
    metrics(i,3) = max(abs(e_alpha(valid)));

    % Beta
    [t_lin, y_lin, ok_lin] = getSignal(pairs{i,3});
    [t_nl,  y_nl,  ok_nl]  = getSignal(pairs{i,4});

    [e_beta, valid] = linearizationError(t_lin, y_lin, t_nl, y_nl);

    metrics(i,2) = sqrt(mean(e_beta(valid).^2));
    metrics(i,4) = max(abs(e_beta(valid)));

end

T_linearization = array2table(metrics, ...
    'VariableNames', {'RMS_alpha_rad','RMS_beta_rad', ...
    'Peak_alpha_rad','Peak_beta_rad'}, ...
    'RowNames', {'mixsyn_hinfsyn','hinfstruct'});

disp(' ');
disp('=== ERRORE DI LINEARIZZAZIONE [rad] ===');
disp(T_linearization);

% Esportazione dei risultati
writetable(T_linearization, ...
    fullfile(outputDir, 'Hinf_linearization_metrics.csv'), ...
    'WriteRowNames', true);

%% ========================================================================
% FUNZIONE LOCALE: errore con interpolazione sul tempo del modello lineare
% ========================================================================

function [e, valid] = linearizationError(t_lin, y_lin, t_nl, y_nl)
    t_lin = t_lin(:);
    y_lin = y_lin(:);
    t_nl  = t_nl(:);
    y_nl  = y_nl(:);
    
    % Rimuove campioni non validi
    idx_lin = isfinite(t_lin) & isfinite(y_lin);
    idx_nl  = isfinite(t_nl)  & isfinite(y_nl);
    
    t_lin = t_lin(idx_lin);
    y_lin = y_lin(idx_lin);
    t_nl  = t_nl(idx_nl);
    y_nl  = y_nl(idx_nl);
    
    % Elimina eventuali tempi duplicati
    [t_lin, ia] = unique(t_lin, 'stable');
    y_lin = y_lin(ia);
    
    [t_nl, ia] = unique(t_nl, 'stable');
    y_nl = y_nl(ia);
    
    % Confronta solo l'intervallo temporale comune
    valid = t_lin >= min(t_nl) & t_lin <= max(t_nl);
    t_common = t_lin(valid);
    
    if isempty(t_common)
        error('I segnali lineare e non lineare non hanno tempi comuni.');
    end
    
    y_nl_interp = interp1(t_nl, y_nl, t_common, 'linear');
    
    e = nan(size(t_lin));
    e(valid) = y_nl_interp - y_lin(valid);

end


%% ========================================================================
% FUNZIONI LOCALI
% ========================================================================

function plotTracking(ctrls,labels,colors,suffix, ...
    angleUnit,tStep,alphaRefDeg,timeVariable)

    % ALPHA
    ax = nexttile;
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    title(ax,['Angolo \alpha - modello ' modelLabel(suffix)]);
    ylabel(ax,'\alpha [deg]');

    h = gobjects(0);
    leg = {};

    for k = 1:numel(ctrls)
        name = ['alpha_' ctrls{k} suffix];
        [t,y,ok] = getSignal(name,timeVariable);

        if ok
            y = angleToDeg(y,angleUnit);
            h(end+1) = plot(ax,t,y,'Color',colors(k,:), ...
                'LineWidth',1.25); %#ok<AGROW>
            leg{end+1} = labels{k}; %#ok<AGROW>
        else
            warning('Variabile non trovata o non utilizzabile: %s',name);
        end
    end

    [t,~,ok] = getSignal(['alpha_' ctrls{1} suffix],timeVariable);

    if ok && numel(t) >= 2
        t0 = min(t);
        tf = max(t);
        ts = min(max(tStep,t0),tf);

        plot(ax,[t0 ts ts tf],[0 0 alphaRefDeg alphaRefDeg], ...
            'k--','LineWidth',1.2);

        h(end+1) = plot(ax,nan,nan,'k--','LineWidth',1.2);
        leg{end+1} = 'r_\alpha';
    end

    if ~isempty(h)
        legend(ax,h,leg,'Location','best');
    end

    % BETA
    ax = nexttile;
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    title(ax,['Angolo \beta - modello ' modelLabel(suffix)]);
    ylabel(ax,'\beta [deg]');
    xlabel(ax,'Tempo [s]');

    h = gobjects(0);
    leg = {};

    for k = 1:numel(ctrls)
        name = ['beta_' ctrls{k} suffix];
        [t,y,ok] = getSignal(name,timeVariable);

        if ok
            y = angleToDeg(y,angleUnit);
            h(end+1) = plot(ax,t,y,'Color',colors(k,:), ...
                'LineWidth',1.25); %#ok<AGROW>
            leg{end+1} = labels{k}; %#ok<AGROW>
        else
            warning('Variabile non trovata o non utilizzabile: %s',name);
        end
    end

    yline(ax,0,'k--','LineWidth',1);
    h(end+1) = plot(ax,nan,nan,'k--');
    leg{end+1} = 'r_\beta = 0';

    legend(ax,h,leg,'Location','best');
end


function plotCommands(ctrls,labels,colors,suffix,timeVariable)
    for channel = 1:2
        ax = nexttile;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        title(ax,sprintf('Comando di controllo u_%d',channel));
        ylabel(ax,sprintf('u_%d [N]',channel));

        if channel == 2
            xlabel(ax,'Tempo [s]');
        end

        h = gobjects(0);
        leg = {};

        for k = 1:numel(ctrls)
            name = ['u_cmd_' ctrls{k} suffix];
            [t,u,ok] = getSignal(name,timeVariable);

            if ~ok
                warning('Variabile non trovata o non utilizzabile: %s',name);
                continue
            end

            if isvector(u)
                % Un vettore contiene un solo canale.
                if channel == 1
                    uc = u(:);
                else
                    continue
                end
            elseif size(u,2) >= channel
                uc = u(:,channel);
            elseif size(u,1) >= channel
                uc = u(channel,:)';
            else
                continue
            end

            h(end+1) = plot(ax,t,uc,'Color',colors(k,:), ...
                'LineWidth',1.15); %#ok<AGROW>
            leg{end+1} = labels{k}; %#ok<AGROW>
        end

        if ~isempty(h)
            legend(ax,h,leg,'Location','best');
        end
    end
end


function plotMeasurementError(ctrls,labels,colors,angleUnit,timeVariable)
    for channel = 1:2
        ax = nexttile;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        if channel == 1
            prefix = 'alpha';
            title(ax,'Errore di misura \alpha');
            ylabel(ax,'e_{\alpha,meas} [deg]');
        else
            prefix = 'beta';
            title(ax,'Errore di misura \beta');
            ylabel(ax,'e_{\beta,meas} [deg]');
            xlabel(ax,'Tempo [s]');
        end

        h = gobjects(0);
        leg = {};

        for k = 1:numel(ctrls)
            c = ctrls{k};

            % Angoli reali: preferenza ai segnali non lineari.
            realNames = {[prefix '_' c '_nl'], ...
                         [prefix '_' c]};

            % Misure: prova suffisso _nl, poi senza suffisso.
            measNames = {[prefix '_meas_' c '_nl'], ...
                         [prefix '_meas_' c]};

            % Compatibilità con il possibile refuso nel nome beta.
            if strcmp(prefix,'beta')
                measNames{end+1} = ['beta_meas_' c 'tuct'];
            end

            [tr,yr,okReal] = getFirstSignal(realNames,timeVariable);
            [tm,ym,okMeas] = getFirstSignal(measNames,timeVariable);

            if ~(okReal && okMeas)
                warning('Misura o angolo reale mancante: %s (%s).',c,prefix);
                continue
            end

            yr = angleToDeg(yr,angleUnit);
            ym = angleToDeg(ym,angleUnit);

            [t,yReal,yMeas] = alignSignals(tr,yr,tm,ym);

            if numel(t) < 2
                continue
            end

            % Errore misura = misura rumorosa - angolo reale.
            err = yMeas-yReal;

            h(end+1) = plot(ax,t,err,'Color',colors(k,:), ...
                'LineWidth',0.9); %#ok<AGROW>
            leg{end+1} = labels{k}; %#ok<AGROW>
        end

        yline(ax,0,'k--');

        if ~isempty(h)
            legend(ax,h,leg,'Location','best');
        end
    end
end


function plotLinearizationError(ctrls,labels,colors,angleUnit,timeVariable)
    for channel = 1:2
        ax = nexttile;
        hold(ax,'on'); grid(ax,'on'); box(ax,'on');

        if channel == 1
            prefix = 'alpha';
            title(ax,'Errore di linearizzazione \alpha');
            ylabel(ax,'\Delta\alpha [deg]');
        else
            prefix = 'beta';
            title(ax,'Errore di linearizzazione \beta');
            ylabel(ax,'\Delta\beta [deg]');
            xlabel(ax,'Tempo [s]');
        end

        h = gobjects(0);
        leg = {};

        for k = 1:numel(ctrls)
            c = ctrls{k};

            [tl,yl,okL] = getSignal([prefix '_' c],timeVariable);
            [tn,yn,okN] = getSignal([prefix '_' c '_nl'],timeVariable);

            if ~(okL && okN)
                warning('Segnali lineare/nonlineare mancanti: %s (%s).',c,prefix);
                continue
            end

            yl = angleToDeg(yl,angleUnit);
            yn = angleToDeg(yn,angleUnit);

            [t,yLinear,yNonlinear] = alignSignals(tl,yl,tn,yn);

            if numel(t) < 2
                continue
            end

            % Differenza tra risposta non lineare e risposta lineare.
            err = yNonlinear-yLinear;

            h(end+1) = plot(ax,t,err,'Color',colors(k,:), ...
                'LineWidth',1.1); %#ok<AGROW>
            leg{end+1} = labels{k}; %#ok<AGROW>
        end

        yline(ax,0,'k--');

        if ~isempty(h)
            legend(ax,h,leg,'Location','best');
        end
    end
end


function T = buildTrackingMetrics(ctrls,angleUnit,tStep,alphaRefDeg,timeVariable)
    Controller = strings(0,1);
    Model = strings(0,1);
    Signal = strings(0,1);
    RMSE = [];
    MAE = [];
    PeakAbsError = [];
    FinalError = [];

    for modelIdx = 1:2
        if modelIdx == 1
            suffix = '';
            model = 'Linear';
        else
            suffix = '_nl';
            model = 'Nonlinear';
        end

        for k = 1:numel(ctrls)
            for channel = 1:2

                if channel == 1
                    prefix = 'alpha';
                    ref = alphaRefDeg;
                else
                    prefix = 'beta';
                    ref = 0;
                end

                [t,y,ok] = getSignal( ...
                    [prefix '_' ctrls{k} suffix],timeVariable);

                if ~ok || numel(t) < 2
                    continue
                end

                y = angleToDeg(y,angleUnit);
                r = ref*ones(size(y));

                if channel == 1
                    r(t < tStep) = 0;
                end

                e = y-r;
                e = e(isfinite(e));

                if isempty(e)
                    continue
                end

                Controller(end+1,1) = string(ctrls{k}); %#ok<AGROW>
                Model(end+1,1) = string(model); %#ok<AGROW>
                Signal(end+1,1) = string(prefix); %#ok<AGROW>

                RMSE(end+1,1) = sqrt(mean(e.^2)); %#ok<AGROW>
                MAE(end+1,1) = mean(abs(e)); %#ok<AGROW>
                PeakAbsError(end+1,1) = max(abs(e)); %#ok<AGROW>
                FinalError(end+1,1) = e(end); %#ok<AGROW>
            end
        end
    end

    T = table(Controller,Model,Signal,RMSE,MAE, ...
        PeakAbsError,FinalError);
end


function [t,y,ok] = getFirstSignal(names,timeVariable)
    t = [];
    y = [];
    ok = false;

    for i = 1:numel(names)
        [t,y,ok] = getSignal(names{i},timeVariable);
        if ok
            return
        end
    end
end



function [t,y,ok] = getSignal(name,timeVariable)
    t = [];
    y = [];
    ok = false;

    % Recupera direttamente la variabile dal workspace base.
    try
        if ischar(name) || (isstring(name) && isscalar(name))

            varName = char(name);

            if ~evalin('base', ...
                    sprintf('exist(''%s'',''var'')', varName))
                return
            end

            s = evalin('base', varName);

        else
            % Il segnale è già stato passato direttamente.
            s = name;
            varName = 'segnale passato direttamente';
        end

    catch ME
        warning('Errore recuperando il segnale: %s', ME.message);
        return
    end

    if isempty(s)
        return
    end

    try
        % ================================================================
        % TIMESERIES
        % ================================================================
        if isa(s,'timeseries')

            t = s.Time(:);
            y = squeeze(s.Data);

            % Un singolo segnale diventa un vettore colonna.
            if isvector(y)
                y = y(:);
            end

        % ================================================================
        % STRUCTURE WITH TIME
        % ================================================================
        elseif isstruct(s) && isfield(s,'time') ...
                && isfield(s,'signals')

            t = s.time(:);
            y = squeeze(s.signals.values);

            if isvector(y)
                y = y(:);
            end

        % ================================================================
        % TIMETABLE
        % ================================================================
        elseif istimetable(s)

            rt = s.Properties.RowTimes;
            t = seconds(rt-rt(1));
            y = table2array(s);

            if isvector(y)
                y = y(:);
            end

        % ================================================================
        % OGGETTI CON PROPRIETA' TIME E DATA
        % ================================================================
        elseif isobject(s) && isprop(s,'Time') ...
                && isprop(s,'Data')

            t = s.Time(:);
            y = squeeze(s.Data);

            if isvector(y)
                y = y(:);
            end

        % ================================================================
        % ARRAY NUMERICI
        % ================================================================
        elseif isnumeric(s)

            y = s;

            if isvector(y)
                y = y(:);
            end

            % Cerca un vettore temporale nel workspace.
            [t,okTime] = findCommonTime(size(y,1),timeVariable);

            if ~okTime
                if isvector(y)
                    t = (0:numel(y)-1)';
                else
                    t = (0:size(y,1)-1)';
                end
            end

        else
            warning('Formato non supportato per %s: %s', ...
                name,class(s));
            return
        end

        % ================================================================
        % CONTROLLO FINALE DELLE DIMENSIONI
        % ================================================================

        % Uniforma le dimensioni del tempo
        t = t(:);

        % Riduce eventuali dimensioni singleton, mantenendo
        % le due colonne dei comandi quando presenti
        y = squeeze(y);

        % Accetta sia segnali scalari sia matrici N x M
        if isvector(y)
            y = y(:);
        end

        % Verifica che il numero di righe coincida con i campioni temporali
        if size(y, 1) == numel(t)
            ok = true;
        elseif size(y, 2) == numel(t)
            % Gestisce anche segnali memorizzati per righe
            y = y.';
            ok = true;
        else
            warning('Dimensioni incompatibili per %s. Tempo: %d campioni, segnale: %s', ...
                name, numel(t), mat2str(size(y)));
            ok = false;
        end

        % Elimina eventuali campioni non validi.
        valid = isfinite(t) & all(isfinite(y), 2);

        t = t(valid);
        y = y(valid, :);

        if numel(t) < 2
            return
        end

        ok = true;

    catch ME
        warning('Errore interpretando %s: %s',name,ME.message);
        t = [];
        y = [];
        ok = false;
    end
end


function [t,ok] = findCommonTime(nSamples,timeVariable)
    t = [];
    ok = false;

    % Prima prova il nome specificato dall'utente.
    candidates = {};

    if ~isempty(timeVariable)
        candidates{end+1} = timeVariable;
    end

    candidates = [candidates, {'tout','t','time'}];

    for i = 1:numel(candidates)
        name = candidates{i};

        try
            if evalin('base',sprintf('exist(''%s'',''var'')',name))
                candidate = evalin('base',name);

                if isnumeric(candidate) && isvector(candidate) ...
                        && numel(candidate) == nSamples

                    t = candidate(:);

                    if all(isfinite(t)) && all(diff(t) >= 0)
                        ok = true;
                        return
                    end
                end
            end
        catch
            % Prova il nome successivo.
        end
    end
end


function tf = isTimeVector(candidate,timeVariable)
    tf = false;

    [t,ok] = findCommonTime(numel(candidate),timeVariable);

    if ok
        % Controllo indicativo: evita di scambiare una colonna qualsiasi
        % per il tempo quando esiste già un vettore temporale comune.
        tf = isequal(t(:),candidate(:));
    end
end


function [t,aAligned,bAligned] = alignSignals(ta,a,tb,b)
    ta = ta(:);
    tb = tb(:);
    a = squeeze(a);
    b = squeeze(b);
    a = a(:);
    b = b(:);

    n = min(numel(ta),numel(a));
    ta = ta(1:n);
    a = a(1:n);

    n = min(numel(tb),numel(b));
    tb = tb(1:n);
    b = b(1:n);

    [ta,ia] = unique(ta,'stable');
    a = a(ia);

    [tb,ib] = unique(tb,'stable');
    b = b(ib);

    t0 = max(min(ta),min(tb));
    tf = min(max(ta),max(tb));

    if tf <= t0
        t = [];
        aAligned = [];
        bAligned = [];
        return
    end

    t = ta(ta >= t0 & ta <= tf);

    if numel(t) < 2
        nGrid = min(1000,max(numel(ta),numel(tb)));
        t = linspace(t0,tf,nGrid)';
    end

    aAligned = interp1(ta,a,t,'linear');
    bAligned = interp1(tb,b,t,'linear');
end


function y = angleToDeg(y,unit)
    if strcmpi(unit,'rad')
        y = y*(180/pi);
    elseif ~strcmpi(unit,'deg')
        error('angleUnit deve essere ''rad'' oppure ''deg''.');
    end
end


function saveFigure(fig,outDir,fileName)
    fullName = fullfile(outDir,fileName);
    exportgraphics(fig,fullName,'Resolution',200);
    fprintf('Salvata figura: %s\n',fullName);
end


function s = onOff(tf)
    if tf
        s = 'on';
    else
        s = 'off';
    end
end


function s = modelLabel(suffix)
    if isempty(suffix)
        s = 'lineare';
    else
        s = 'non lineare';
    end
end


function plotMaxSigma(sys1, sys2, w, label1, label2)
% Confronta i massimi valori singolari di due sistemi MIMO.
% I risultati sono espressi in dB.
    
    sig1 = maxSigmaDb(sys1, w);
    sig2 = maxSigmaDb(sys2, w);
    
    semilogx(w, sig1, 'b-', 'LineWidth', 1.5);
    hold on;
    semilogx(w, sig2, 'r--', 'LineWidth', 1.5);
    hold off;
    
    grid on;
    xlim([w(1) w(end)]);
    ylim([-80 50]);
    
    ylabel('Modulo [dB]');
    legend(label1, label2, 'Location','best');
end


function ydb = maxSigmaDb(sys, w)
% Calcola il massimo valore singolare di sys a ogni frequenza.
    
    H = freqresp(sys, w);
    ydb = zeros(1, numel(w));

    for k = 1:numel(w)
        Hk = H(:,:,k);
        sv = svd(Hk);
        ydb(k) = 20*log10(max(sv));
    end
end
