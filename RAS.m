%% NOTE TEORICHE - REGIONE DI STABILITA' ASINTOTICA
% La RAS viene stimata numericamente su un orizzonte finito. Un punto e'
% classificato nella regione se, senza violare i limiti globali, rimane
% entro le tolleranze di stato durante la coda finale della simulazione.
% La mappa e' quindi una stima numerica conservativa della regione di
% attrazione, dipendente da tFinal e dalle tolleranze dichiarate.
%
% =========================================================================
% RAS numerica vettorializzata - closed-loop CONTINUO
%
% V4:
%   - plant non lineare: RK4 vettorializzato;
%   - ctrl lineari: propagazione ESATTA via matrice esponenziale
%     per ingresso sensor-ZOH costante nel passo;
%   - motori EMAX: RK4 continuo, pilotati dal comando del ctrl
%     valutato a t, t+h/2 e t+h;
%   - saturazione a ogni stadio del motore;
%   - sensori VN-100: ZOH;
%   - transport delay: buffer campionato corretto ai sample time.
%
% La propagazione esatta dei ctrl evita l'instabilità numerica
% dell'RK4 esplicito sul controllore mu di ordine elevato.
%
% MATLAB R2025a / R2026a.
% =========================================================================
close all;
clc;

%% ========================================================================
% 0. CONFIGURAZIONE
% ========================================================================
opt.alphaRangeDeg = [-180 180];
opt.betaRangeDeg  = [-180 180];


opt.nAlpha = 211;
opt.nBeta  = 321;
opt.tFinal = 30;             % [s] orizzonte della stima numerica RAS
% 1 ms divide esattamente il delay nominale di 15 ms.
% Il passo viene scelto sufficientemente piccolo per la propagazione
% RK4 del plant e dei motori e per la gestione del transport delay.
opt.dt = 0.001;             % [s]
opt.tailFraction = 0.10;
opt.angleToleranceDeg = 0.25;
opt.rateToleranceDeg  = 0.50;

opt.maxAlphaDeg = 85;
opt.maxBetaDeg  = 170;
opt.earlyReject = true;
opt.saveMAT = true;
opt.saveFigures = false;
opt.ctrlNames = { ...
    'LQG1', ...
   %% 'LQG', ...
   %% 'LQGI', ...
   %% 'mixsyn', ...
   %% 'hinfsyn', ...
   %% 'PID+comp', ...
   %% 'mu-synthesis'
};

%% ========================================================================
% 1. INIZIALIZZAZIONE DEL PROGETTO
% ========================================================================
if ~exist('p','var') || ~exist('act','var') || ...
        ~exist('sensor','var') || ~exist('alpha_0','var') || ...
        ~exist('beta_0','var') || ~exist('u0_nom','var') || ~exist('x0','var')
    dataset_elicottero; % Chiama il tuo script
end
% RAS deterministica
sensor.noiseEnable = 0;
act.noiseEnable = 0;


%% ========================================================================
% 2. CARICAMENTO DEI CONTROLLORI CONTINUI
% ========================================================================
ctrls = loadctrls(opt.ctrlNames);
fprintf('\n============================================================\n');
fprintf('RAS FINAL V5 - CONTROLLORI CARICATI\n');
fprintf('============================================================\n');
for k = 1:numel(ctrls)
    fprintf('%-16s ordine = %d, ingressi = %d, uscite = %d\n', ...
        ctrls(k).name, ...
        size(ctrls(k).A,1), ...
        size(ctrls(k).B,2), ...
        size(ctrls(k).C,1));
end
% Transizioni esatte del ctrl per ingresso costante nel passo.
% I sensori sono ZOH, quindi uk è effettivamente costante tra due update.
for k = 1:numel(ctrls)
    nk = size(ctrls(k).A,1);
    if nk > 0
        [ctrls(k).Phi,ctrls(k).Gamma] = ...
            exactZOH(ctrls(k).A,ctrls(k).B,opt.dt);
        [ctrls(k).PhiHalf,ctrls(k).GammaHalf] = ...
            exactZOH(ctrls(k).A,ctrls(k).B,opt.dt/2);
    else
        ctrls(k).Phi = zeros(0);
        ctrls(k).Gamma = zeros(0,size(ctrls(k).B,2));
        ctrls(k).PhiHalf = zeros(0);
        ctrls(k).GammaHalf = zeros(0,size(ctrls(k).B,2));
    end
end

%% ========================================================================
% 3. SENSORI E JACOBIANO (Allineato a dataset_elicottero.m)
% ========================================================================
% Matrice della varianza dei sensori
V_sensor = diag([sensor.acc.var, sensor.mag.var, sensor.mag.var]);
Vinv = diag(1./diag(V_sensor));

% Valore di equilibrio dei sensori (il tuo y0)
y0_sensor = [
    -p.g * sin(alpha_0);
     cos(alpha_0) * cos(beta_0);
    -sin(beta_0)
];

% Jacobiano dei sensori calcolato nel punto di equilibrio
Jsensor = [
    -p.g*cos(alpha_0),            0;
    -sin(alpha_0)*cos(beta_0),     -cos(alpha_0)*sin(beta_0);
     0,                          -cos(beta_0)
];

Hy = (Jsensor.'*Vinv*Jsensor) \ (Jsensor.'*Vinv);

%% ========================================================================
% 4. DINAMICA CONTINUA DEI MOTORI E TRANSPORT DELAY
% ========================================================================
Am1 = [
    0, 1;
   -act.wn1^2, -2*act.zeta1*act.wn1
];
Bm1 = [0; act.wn1^2];
Cm1 = [1 0];
Am2 = [
    0, 1;
   -act.wn2^2, -2*act.zeta2*act.wn2
];
Bm2 = [0; act.wn2^2];
Cm2 = [1 0];
delaySteps1 = round(act.td1/opt.dt);
delaySteps2 = round(act.td2/opt.dt);
delayError1 = abs(delaySteps1*opt.dt-act.td1);
delayError2 = abs(delaySteps2*opt.dt-act.td2);
if delayError1 > 1e-12 || delayError2 > 1e-12
    warning(['opt.dt non divide esattamente il transport delay. ', ...
             'Il solo ritardo viene quantizzato alla griglia temporale.']);
end
delaySteps1 = max(delaySteps1,1);
delaySteps2 = max(delaySteps2,1);
fprintf('\nTransport delay 1: %d campioni = %.6f s\n', ...
    delaySteps1,delaySteps1*opt.dt);
fprintf('Transport delay 2: %d campioni = %.6f s\n', ...
    delaySteps2,delaySteps2*opt.dt);

%% ========================================================================
% 5. GRIGLIA REGOLARE DELLE CONDIZIONI INIZIALI
% ========================================================================
[deltaAlpha0,deltaBeta0,alphaAxis,betaAxis,nTraj] = ...
    buildInitialGrid(opt);
fprintf('\n============================================================\n');
fprintf('GRIGLIA RAS FINAL V5\n');
fprintf('============================================================\n');
fprintf('alpha: %d punti\n',opt.nAlpha);
fprintf('beta : %d punti\n',opt.nBeta);
fprintf('totale per controllore: %d traiettorie\n',nTraj);
fprintf('Integrazione closed-loop vettorializzata continua.\n');

%% ========================================================================
% 6. SOGLIE
% ========================================================================
angleTol = deg2rad(opt.angleToleranceDeg);
rateTol  = deg2rad(opt.rateToleranceDeg);
maxAlpha = deg2rad(opt.maxAlphaDeg);
maxBeta  = deg2rad(opt.maxBetaDeg);
tailStartTime = opt.tFinal*(1-opt.tailFraction);
nIter = round(opt.tFinal/opt.dt);

%% ========================================================================
% 7. CALCOLO RAS
% ========================================================================
results = struct([]);
for ic = 1:numel(ctrls)
    Ctl = ctrls(ic);
    fprintf('\n============================================================\n');
    fprintf('RAS FINAL V5: %s\n',upper(Ctl.name));
    fprintf('============================================================\n');
    tic;
    % Stato fisico: [alpha; alpha_dot; beta; beta_dot]
    X = repmat(x0,1,nTraj);
    X(1,:) = alpha_0 + deltaAlpha0;
    X(3,:) = beta_0  + deltaBeta0;
    % Stato continuo del controllore
    nk = size(Ctl.A,1);
    Xk = zeros(nk,nTraj);
    % Stati continui dei motori
    Xm1 = zeros(2,nTraj);
    Xm2 = zeros(2,nTraj);
    % Delay buffer dell'uscita dei motori
    % +1 perché il plant al tempo t_n deve leggere y(t_n-td), non
    % il campione già avanzato di un passo.
    delayBuffer1 = zeros(delaySteps1+1,nTraj);
    delayBuffer2 = zeros(delaySteps2+1,nTraj);
    delayPtr1 = 1;
    delayPtr2 = 1;
    % Sensori campionati/ZOH
    alphaSensor = X(1,:);
    betaSensor  = X(3,:);
    accEvery = max(1,round(sensor.acc.Ts/opt.dt));
    magEvery = max(1,round(sensor.mag.Ts/opt.dt));
    alive = true(1,nTraj);
    discarded = false(1,nTraj);
    tailErrorMax = zeros(4,nTraj);
    tailStartErrorNorm = nan(1,nTraj);
    maxAbsAlpha = abs(X(1,:));
    maxAbsBeta  = abs(X(3,:));
    maxAbsCommand1 = zeros(1,nTraj);
    maxAbsCommand2 = zeros(1,nTraj);
    for istep = 1:nIter
        if ~any(alive)
            break
        end
        traj = find(alive);
        % -------------------------------------------------------------
        % Sensori deterministici con ZOH
        % -------------------------------------------------------------
        if mod(istep-1,accEvery) == 0
            alphaSensor(traj) = X(1,traj);
        end
        if mod(istep-1,magEvery) == 0
            betaSensor(traj) = X(3,traj);
        end
        y = [
            -p.g * sin(alphaSensor(traj));
             cos(alphaSensor(traj)) .* cos(betaSensor(traj));
            -sin(betaSensor(traj))
        ];
        deltaY = y - y0_sensor;
% -------------------------------------------------------------
        % Ingresso al ctrl, mantenuto costante durante il passo RK4
        % -------------------------------------------------------------
        if size(Ctl.B, 2) == 5
            % LQG 2-DOF o con Integratore: ingressi [r_alpha; r_beta; y_acc; mx; my]
            % Il riferimento r è nullo perché valutiamo la stabilità attorno all'equilibrio
            uk = [zeros(2,numel(traj)); deltaY];
            
        elseif size(Ctl.B, 2) == 3
            % LQG 1-DOF Regolatore: ingressi [y_acc; mx; my]
            uk = deltaY;
            
        else
            % Controllori H-inf (Mixsyn, Hinfstruct): ingressi [e_alpha; e_beta]
            % dove e = r - y_stimato. Essendo r=0, e = -y_stimato
            deltaAngleEstimated = Hy*deltaY;
            uk = -deltaAngleEstimated;
        end
        % -------------------------------------------------------------
        % Forze DELAYED viste dal plant durante il passo corrente.
        % Il transport delay rende questa uscita già dipendente dal passato.
        % -------------------------------------------------------------
        deltaF1Delayed = delayBuffer1(delayPtr1,traj);
        deltaF2Delayed = delayBuffer2(delayPtr2,traj);
        Fdelayed = [
            u0_nom(1) + deltaF1Delayed;
            u0_nom(2) + deltaF2Delayed
        ];
        % -------------------------------------------------------------
        % RK4 DEL PLANT NON LINEARE
        % -------------------------------------------------------------
        X(:,traj) = rk4Step(X(:,traj),Fdelayed,p,opt.dt);
        % -------------------------------------------------------------
        % ctrl: PROPAGAZIONE ESATTA; MOTORI: RK4 CONTINUO
        %
        % uk è costante nel passo per effetto degli ZOH sensore.
        % Calcoliamo esattamente lo stato del ctrl a:
        %   t, t+h/2, t+h
        % e usiamo i relativi comandi per gli stadi RK4 dei motori.
        % -------------------------------------------------------------
        Xk0 = Xk(:,traj);
        if nk > 0
            XkHalf = ...
                Ctl.PhiHalf*Xk0 + ...
                Ctl.GammaHalf*uk;
            XkEnd = ...
                Ctl.Phi*Xk0 + ...
                Ctl.Gamma*uk;
        else
            XkHalf = zeros(0,numel(traj));
            XkEnd  = zeros(0,numel(traj));
        end
        [cmd1_1,cmd2_1] = ctrlCommand(Xk0,uk,Ctl,act);
        [cmd1_2,cmd2_2] = ctrlCommand(XkHalf,uk,Ctl,act);
        [cmd1_4,cmd2_4] = ctrlCommand(XkEnd,uk,Ctl,act);
        % Motore 1
        m1_1 = Am1*Xm1(:,traj) + Bm1*cmd1_1;
        m1_2 = Am1*(Xm1(:,traj)+0.5*opt.dt*m1_1) + Bm1*cmd1_2;
        m1_3 = Am1*(Xm1(:,traj)+0.5*opt.dt*m1_2) + Bm1*cmd1_2;
        m1_4 = Am1*(Xm1(:,traj)+opt.dt*m1_3)     + Bm1*cmd1_4;
        Xm1(:,traj) = Xm1(:,traj) + ...
            (opt.dt/6)*(m1_1 + 2*m1_2 + 2*m1_3 + m1_4);
        % Motore 2
        m2_1 = Am2*Xm2(:,traj) + Bm2*cmd2_1;
        m2_2 = Am2*(Xm2(:,traj)+0.5*opt.dt*m2_1) + Bm2*cmd2_2;
        m2_3 = Am2*(Xm2(:,traj)+0.5*opt.dt*m2_2) + Bm2*cmd2_2;
        m2_4 = Am2*(Xm2(:,traj)+opt.dt*m2_3)     + Bm2*cmd2_4;
        Xm2(:,traj) = Xm2(:,traj) + ...
            (opt.dt/6)*(m2_1 + 2*m2_2 + 2*m2_3 + m2_4);
        if nk > 0
            Xk(:,traj) = XkEnd;
        end
        maxAbsCommand1(traj) = max( ...
            maxAbsCommand1(traj), ...
            max([abs(cmd1_1); abs(cmd1_2); abs(cmd1_4)],[],1));
        maxAbsCommand2(traj) = max( ...
            maxAbsCommand2(traj), ...
            max([abs(cmd2_1); abs(cmd2_2); abs(cmd2_4)],[],1));
        % -------------------------------------------------------------
        % Nuova uscita dei motori nel delay buffer
        % -------------------------------------------------------------
        motorOut1 = Cm1*Xm1(:,traj);
        motorOut2 = Cm2*Xm2(:,traj);
        delayBuffer1(delayPtr1,traj) = motorOut1;
        delayBuffer2(delayPtr2,traj) = motorOut2;
        delayPtr1 = delayPtr1 + 1;
        if delayPtr1 > (delaySteps1+1)
            delayPtr1 = 1;
        end
        delayPtr2 = delayPtr2 + 1;
        if delayPtr2 > (delaySteps2+1)
            delayPtr2 = 1;
        end
        % -------------------------------------------------------------
        % Diagnostica / early reject
        % -------------------------------------------------------------
        maxAbsAlpha(traj) = max(maxAbsAlpha(traj),abs(X(1,traj)));
        maxAbsBeta(traj)  = max(maxAbsBeta(traj), abs(X(3,traj)));
        finiteNow = ...
            all(isfinite(X(:,traj)),1) & ...
            all(isfinite(Xm1(:,traj)),1) & ...
            all(isfinite(Xm2(:,traj)),1);
        if nk > 0
            finiteNow = finiteNow & all(isfinite(Xk(:,traj)),1);
        end
        badAngle = checkAngleLimits(X(:,traj),maxAlpha,maxBeta);
        badLocal = ~finiteNow | badAngle;
        if any(badLocal)
            badtraj = traj(badLocal);
            discarded(badtraj) = true;
            if opt.earlyReject
                alive(badtraj) = false;
            end
        end
        % -------------------------------------------------------------
        % Criterio di convergenza nell'ultimo 10 %
        % -------------------------------------------------------------
        tNext = istep*opt.dt;
        if tNext >= tailStartTime
            trajTail = find(alive & ~discarded);
            if ~isempty(trajTail)
                err = abs(X(:,trajTail)-x0);

                firstTail = isnan(tailStartErrorNorm(trajTail));
                if any(firstTail)
                    trajFirst = trajTail(firstTail);
                    tailStartErrorNorm(trajFirst) = ...
                        vecnorm(X(:,trajFirst)-x0,2,1);
                end

                tailErrorMax(:,trajTail) = max( ...
                    tailErrorMax(:,trajTail),err);
            end
        end
        if mod(istep,max(1,round(nIter/10))) == 0
            fprintf('  %5.1f %% | attive: %d / %d\n', ...
                100*istep/nIter,nnz(alive),nTraj);
        end
    end

    %% Classificazione finale
    accepted = classifyRAS( ...
        discarded, ...
        maxAbsAlpha, ...
        maxAbsBeta, ...
        tailErrorMax, ...
        angleTol, ...
        rateTol, ...
        maxAlpha, ...
        maxBeta);

    finalErrorNorm = vecnorm(X-x0,2,1);
    validTrajectory = ...
        ~discarded & ...
        maxAbsAlpha < maxAlpha & ...
        maxAbsBeta  < maxBeta;
    
    slowConverging = ...
        validTrajectory & ...
        ~accepted & ...
        isfinite(tailStartErrorNorm) & ...
        (finalErrorNorm < tailStartErrorNorm);

    attractionMap = reshape(accepted,opt.nBeta,opt.nAlpha);
    elapsed = toc;
    fprintf('\n%s completato in %.2f s\n',Ctl.name,elapsed);
    fprintf('Punti stabili: %d / %d (%.1f %%)\n', ...
        nnz(accepted),nTraj,100*nnz(accepted)/nTraj);
    results(ic).name = Ctl.name;
    results(ic).family = Ctl.family;
    results(ic).accepted = accepted;
    results(ic).attractionMap = attractionMap;
    results(ic).tailErrorMax = tailErrorMax;
    results(ic).tailStartErrorNorm = tailStartErrorNorm;
    results(ic).finalErrorNorm = finalErrorNorm;
    results(ic).slowConverging = slowConverging;
    results(ic).maxAbsAlpha = maxAbsAlpha;
    results(ic).maxAbsBeta = maxAbsBeta;
    results(ic).maxAbsCommand1 = maxAbsCommand1;
    results(ic).maxAbsCommand2 = maxAbsCommand2;
    results(ic).elapsed = elapsed;
    results(ic).domainFullyStable = all(accepted);
    results(ic).domainFullyUnstable = ~any(accepted);
end

%% ========================================================================
% 8. PLOT RAS - VERA REGIONE SU GRIGLIA REGOLARE
% ========================================================================
for ic = 1:numel(results)
    Mplot = logical(results(ic).attractionMap);
    figure( ...
        'Name',['RAS FINAL V5 - ',results(ic).name], ...
        'Color','w');
    imagesc( ...
        alphaAxis, ...
        betaAxis, ...
        double(Mplot));
    axis xy
    axis equal
    hold on
    grid on
    box on
    if any(Mplot(:)) && ~all(Mplot(:))
        contour( ...
            alphaAxis, ...
            betaAxis, ...
            double(Mplot), ...
            [0.5 0.5], ...
            'k', ...
            'LineWidth',1.8);
        domainNote = '';
    elseif all(Mplot(:))
        domainNote = '  [RAS oltre il dominio esplorato]';
    else
        domainNote = '  [nessun punto convergente entro il criterio]';
    end
    plot(0,0,'k+','MarkerSize',10,'LineWidth',2);
    xlabel('$\Delta\alpha(0)\;[^{\circ}]$','Interpreter','latex');
    ylabel('$\Delta\beta(0)\;[^{\circ}]$','Interpreter','latex');
    title([ ...
        'Stima numerica della regione di attrazione - ', ...
        results(ic).name, ...
        domainNote], ...
        'Interpreter','latex');
    xlim(opt.alphaRangeDeg);
    ylim(opt.betaRangeDeg);
    cb = colorbar;
    cb.Ticks = [0 1];
    cb.TickLabels = {'Fuori criterio a T_f','Convergente entro T_f'};
    if opt.saveFigures
        exportgraphics( ...
            gcf, ...
            ['RAS_FINAL_V5_',sanitizeName(results(ic).name),'.pdf'], ...
            'ContentType','vector');
    end
end

%% ========================================================================
% 9. FIGURA COMPARATIVA DEI CONFINI
% ========================================================================
figure( ...
    'Name','Confronto confini RAS FINAL V5', ...
    'Color','w');
hold on
grid on
box on
axis equal
hLegend = gobjects(0);

%% Colore diverso per ciascun controllore
colors = lines(numel(results));
for ic = 1:numel(results)
    Mplot = logical(results(ic).attractionMap);
    if any(Mplot(:)) && ~all(Mplot(:))
        [~,h] = contour( ...
            alphaAxis, ...
            betaAxis, ...
            double(Mplot), ...
            [0.5 0.5], ...
            'LineWidth',1.8, ...
            'DisplayName',results(ic).name);
        h.LineColor = colors(ic,:);
        hLegend(end+1) = h; %#ok<SAGROW>
    end
end
plot( ...
    0,0, ...
    'k+', ...
    'MarkerSize',10, ...
    'LineWidth',2, ...
    'HandleVisibility','off');
xlabel('$\Delta\alpha(0)\;[^{\circ}]$','Interpreter','latex');
ylabel('$\Delta\beta(0)\;[^{\circ}]$','Interpreter','latex');
title('Confronto delle stime numeriche delle regioni di attrazione','Interpreter','latex');
xlim(opt.alphaRangeDeg);
ylim(opt.betaRangeDeg);
if ~isempty(hLegend)
    legend( ...
        hLegend, ...
        'Location','bestoutside');
end

%% ========================================================================
% 10. RIEPILOGO
% ========================================================================
fprintf('\n============================================================\n');
fprintf('RIEPILOGO RAS FINAL V5\n');
fprintf('============================================================\n');
for ic = 1:numel(results)
    if results(ic).domainFullyStable
        extra = ' - bordo NON ancora raggiunto';
    elseif results(ic).domainFullyUnstable
        extra = ' - nessun punto convergente entro il criterio';
    else
        extra = ' - confine interno al dominio';
    end
    fprintf('%-16s : %5d / %5d nella RAS numerica (%.1f%%), %5d convergenti lenti%s\n', ...
        results(ic).name, ...
        nnz(results(ic).accepted), ...
        nTraj, ...
        100*nnz(results(ic).accepted)/nTraj, ...
        nnz(results(ic).slowConverging), ...
        extra);
end

%% ========================================================================
% 11. METRICHE QUANTITATIVE DELLA RAS
% ========================================================================
dAlpha = abs(alphaAxis(2)-alphaAxis(1));
dBeta  = abs(betaAxis(2)-betaAxis(1));
cellAreaDeg2 = dAlpha*dBeta;
ctrl = strings(numel(results),1);
StablePoints = zeros(numel(results),1);
StablePercent = zeros(numel(results),1);
SlowConvergingPoints = zeros(numel(results),1);
AreaDeg2 = zeros(numel(results),1);
DeltaAlphaMinDeg = nan(numel(results),1);
DeltaAlphaMaxDeg = nan(numel(results),1);
DeltaBetaMinDeg = nan(numel(results),1);
DeltaBetaMaxDeg = nan(numel(results),1);
TouchesSearchBoundary = false(numel(results),1);
for ic = 1:numel(results)
    M = logical(results(ic).attractionMap);
    ctrl(ic) = string(results(ic).name);
    StablePoints(ic) = nnz(M);
    StablePercent(ic) = 100*nnz(M)/numel(M);
    SlowConvergingPoints(ic) = nnz(results(ic).slowConverging);
    % Area numerica di cella. Per una griglia sufficientemente fitta è una
    % stima semplice, trasparente e riproducibile.
    AreaDeg2(ic) = nnz(M)*cellAreaDeg2;
    if any(M(:))
        [ib,ia] = find(M);
        DeltaAlphaMinDeg(ic) = min(alphaAxis(ia));
        DeltaAlphaMaxDeg(ic) = max(alphaAxis(ia));
        DeltaBetaMinDeg(ic)  = min(betaAxis(ib));
        DeltaBetaMaxDeg(ic)  = max(betaAxis(ib));
        edgeMask = [ ...
            M(1,:), ...
            M(end,:), ...
            M(:,1).', ...
            M(:,end).' ];
        
        TouchesSearchBoundary(ic) = any(edgeMask);
    end
end
RASmetrics = table( ...
    ctrl, ...
    StablePoints, ...
    StablePercent, ...
    SlowConvergingPoints, ...
    AreaDeg2, ...
    DeltaAlphaMinDeg, ...
    DeltaAlphaMaxDeg, ...
    DeltaBetaMinDeg, ...
    DeltaBetaMaxDeg, ...
    TouchesSearchBoundary);
disp(' ');
disp('============================================================');
disp('METRICHE RAS FINALI');
disp('============================================================');
disp(RASmetrics);
writetable(RASmetrics,'RAS_metrics_final_v5.csv');

%% ========================================================================
% 11. SALVATAGGIO
% ========================================================================
if opt.saveMAT
    save( ...
        'RAS_vectorized_results_final_v5.mat', ...
        'results', ...
        'opt', ...
        'alphaAxis', ...
        'betaAxis', ...
        'RASmetrics');
    fprintf('\nSalvato: RAS_vectorized_results_final_v5.mat\n');
end
disp('RAS FINAL V5 completata.');

%% ========================================================================
% FUNZIONI LOCALI
% ========================================================================
function ctrls = loadctrls(reqNames)

    ctrls = struct( ...
        'name',{}, ...
        'family',{}, ...
        'A',{}, ...
        'B',{}, ...
        'C',{}, ...
        'D',{});

    %% ================================================================
    % LQG
    % ================================================================

    if isfile('LQG_Controllers.mat')
        SL = load('LQG_Controllers.mat');
    else
        SL = struct;
    end
    
    % LQG 1-DOF
    if ismember('LQG1', reqNames)
        K = fetchSystem('K_lqg_reg', SL);
        ctrls(end+1) = ...
            makectrl('LQG 1-DOF', 'LQG', K);
    end
    
    % LQG 2-DOF senza integratore
    if ismember('LQG', reqNames)
        K = fetchSystem('K_lqg_2dof', SL);
        ctrls(end+1) = ...
            makectrl('LQG 2-DOF', 'LQG', K);
    end
    
    % LQG 2-DOF con integratore
    if ismember('LQGI', reqNames)
        K = fetchSystem('K_lqg_int', SL);
        ctrls(end+1) = ...
            makectrl('LQG con Int', 'LQG', K);
    end

    %% ================================================================
    % H-INFINITY
    % ================================================================

    if isfile('HINF_controllers.mat')
        SH = load('HINF_controllers.mat');
    else
        SH = struct;
    end

    % Mixsyn
    if ismember('mixsyn', reqNames)
        K = fetchSystem('K_mix', SH);
        ctrls(end+1) = ...
            makectrl('mixsyn', 'HINF', K);
    end

    % Hinfstruct
    if ismember('hinfsyn', reqNames)
        K = fetchSystem('K_hinfsyn', SH);
        ctrls(end+1) = ...
            makectrl('hinfsyn', 'HINF', K);
    end

    % PID + compensatore
    if ismember('PID+comp', reqNames)
        K = fetchSystem('K_pidcomp', SH);
        ctrls(end+1) = ...
            makectrl('PID+comp', 'HINF', K);
    end


    %% ================================================================
    % MU-SYNTHESIS
    % ================================================================

    if ismember('mu-synthesis', reqNames)

        if isfile('MU_controller.mat')
            SM = load('MU_controller.mat');
        else
            SM = struct;
        end

        K = fetchSystem('K_mu', SM);

        ctrls(end+1) = ...
            makectrl('mu-synthesis', 'HINF', K);
    end


    %% ================================================================
    % H2
    % ================================================================

    if ismember('H2', reqNames)

        if isfile('H2_controller.mat')
            S2 = load('H2_controller.mat');
        else
            S2 = struct;
        end

        K = fetchSystem('K_h2', S2);

        ctrls(end+1) = ...
            makectrl('H2', 'HINF', K);
    end


    %% ================================================================
    % CONTROLLO FINALE
    % ================================================================

    if isempty(ctrls)
        error( ...
            ['Nessun controllore richiesto trovato. ', ...
             'Controllare opt.ctrlNames e i file .mat.']);
    end

end

function [dA,dB,aGrid,bGrid,N] = buildInitialGrid(opt)

    aGrid = linspace(opt.alphaRangeDeg(1), ...
                     opt.alphaRangeDeg(2),opt.nAlpha);

    bGrid = linspace(opt.betaRangeDeg(1), ...
                     opt.betaRangeDeg(2),opt.nBeta);

    [AA,BB] = meshgrid(aGrid,bGrid);

    dA = deg2rad(AA(:)).';
    dB = deg2rad(BB(:)).';

    N = numel(dA);
end

function K = fetchSystem(varName, S)

    % Prima cerca nel file .mat
    if isfield(S, varName)
        K = S.(varName);
        return
    end

    % Se non c'è nel file, cerca nel Workspace
    if evalin('base', sprintf("exist('%s','var')", varName))
        K = evalin('base', varName);
        return
    end

    error( ...
        'Controllore %s non trovato nel file .mat né nel Workspace.', ...
        varName);

end

function Ctl = makectrl(name, family, K)
    K_ss = ss(K);
    [A, B, C, D] = ssdata(K_ss);
    Ctl.name = name;
    Ctl.family = family;
    Ctl.A = double(A);
    Ctl.B = double(B);
    Ctl.C = double(C);
    Ctl.D = double(D);
end

function [Phi,Gamma] = exactZOH(A,B,h)
    % Exact state transition for xdot=A*x+B*u with u constant over h.
    n = size(A,1);
    m = size(B,2);
    if n == 0
        Phi = zeros(0);
        Gamma = zeros(0,m);
        return
    end
    E = expm([A, B; zeros(m,n+m)]*h);
    Phi = E(1:n,1:n);
    Gamma = E(1:n,n+(1:m));
end

function dX = helicopterDerivativeVectorized(X,F,p)
    alpha    = X(1,:);
    alphaDot = X(2,:);
    beta     = X(3,:);
    betaDot  = X(4,:);
    F1 = F(1,:);
    F2 = F(2,:);
    Jbeta = ...
        p.J_y.*sin(alpha).^2 + ...
        (p.J_z + p.m*p.l^2).*cos(alpha).^2 + ...
        p.I_b;
    tauAlpha = ...
        p.l.*( ...
            cos(beta).*F1 + ...
            p.eps_p.*sin(beta).*F2) ...
        - p.c_alpha.*alphaDot ...
        - p.m*p.g*p.l.*sin(alpha);
    tauBeta = ...
        p.l.*( ...
            p.eps_y.*sin(alpha).*F1 + ...
            cos(alpha).*F2) ...
        - p.c_beta.*betaDot;
    alphaDDot = tauAlpha ./ p.J_alpha;
    betaDDot  = tauBeta  ./ Jbeta;
    dX = [
        alphaDot;
        alphaDDot;
        betaDot;
        betaDDot
    ];
end

function [cmd1,cmd2] = ctrlCommand(Xk,uk,Ctl,act)
    if isempty(Xk)
        raw = Ctl.D*uk;
    else
        raw = Ctl.C*Xk + Ctl.D*uk;
    end
    cmd1 = min(max( ...
        raw(1,:), ...
        act.deltaF1_min), ...
        act.deltaF1_max);
    cmd2 = min(max( ...
        raw(2,:), ...
        act.deltaF2_min), ...
        act.deltaF2_max);
end

function out = sanitizeName(in)
    out = regexprep(in,'[^a-zA-Z0-9_-]','_');
end
function bad = checkAngleLimits(X,maxA,maxB)

    bad = abs(X(1,:)) >= maxA | ...
          abs(X(3,:)) >= maxB;

end

function accepted = classifyRAS( ...
    discarded,maxAbsAlpha,maxAbsBeta,tailErrorMax,...
    angleTol,rateTol,maxAlpha,maxBeta)

    validAngles = ...
        maxAbsAlpha < maxAlpha & ...
        maxAbsBeta  < maxBeta;

    settled = all( ...
        tailErrorMax([1 3],:) < angleTol & ...
        tailErrorMax([2 4],:) < rateTol, ...
        1);

    accepted = ~discarded & validAngles & settled;

end
function Xnext = rk4Step(X,F,p,h)

    k1 = helicopterDerivativeVectorized(X,F,p);
    k2 = helicopterDerivativeVectorized(X + 0.5*h*k1,F,p);
    k3 = helicopterDerivativeVectorized(X + 0.5*h*k2,F,p);
    k4 = helicopterDerivativeVectorized(X + h*k3,F,p);

    Xnext = X + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end