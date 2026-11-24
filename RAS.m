% =========================================================================
% RAS numerica (non lineare + lineare di controllo) - Elicottero 2DoF
% Controllori: LQG 1-DOF / LQG 2-DOF / LQGI / mixsyn / hinfsyn / PID+comp /
%              mu-synthesis
close all;
clc;

%% ========================================================================
% 0. CONFIGURAZIONE
% ========================================================================
opt.quick          = false;   % true = prova rapida (griglie 9x9, Tf corto)
opt.dt             = 0.002;   % passo [s]. Il delay (0.02 s) deve essere multiplo
opt.ctrlNames      = {'LQG1','LQG2','LQGI','mixsyn','hinfsyn','PID+comp','mu-synthesis'};
opt.lqgRows        = 1:2;     % righe di uscita di K usate come comando [F1;F2]
opt.skipDuplicates = true;    % non ri-simula controllori identici
opt.muReduceOrder  = 0;       % >0: tenta balred del controllore mu a questo ordine
opt.selfTest       = true;    % confronto lineare vs non lineare su 0.2 deg
opt.earlyStop      = true;    % arresta ogni traiettoria appena converge
opt.f1NonNegative  = false;   % true: limita anche la forza ASSOLUTA F1 >= 0
opt.saveMAT        = true;

% --- Criterio di convergenza (stato completo) ---
tol.convAngleDeg   = 1e-3;    % |d_alpha|,|d_beta| < 1e-3 deg
tol.convRateDeg    = 1e-3;    % |velocita'| < 1e-3 deg/s
tol.relBlock       = 1e-3;    % attuatori, ritardo, controllore <= relBlock * proprio picco
tol.calmWindow     = 5;       % [s] le condizioni sopra devono valere in modo continuo
tol.checkEvery     = 0.5;     % [s] intervallo di controllo della convergenza
% --- Classe 2 (indeterminata): decadimento della norma nella coda ---
tol.decayFactor    = 0.999;   % max norma 2a meta coda < decayFactor * 1a meta
tol.tailFraction   = 0.25;    % ultima frazione di Tf usata come "coda"
% --- Limiti fisici (fuori da questi la traiettoria e' scartata) ---
tol.maxAlphaDeg    = 85;
tol.maxBetaDeg     = 170;

% --- Verifica finale ---
verify.enable      = true;
verify.factor      = 3;       % indeterminati ricalcolati con Tf = factor * Tf campagna
verify.nSample     = 20;      % punti campionati per gruppo

% --- Campagne ---
camp = struct('name',{},'alphaRangeDeg',{},'betaRangeDeg',{}, ...
              'nAlpha',{},'nBeta',{},'refine',{},'tFinal',{},'modes',{},'satCases',{});

camp(1).name          = 'LOCALE';
camp(1).alphaRangeDeg = [-30 30];
camp(1).betaRangeDeg  = [-30 30];
camp(1).nAlpha        = 41;
camp(1).nBeta         = 41;
camp(1).refine        = 4;      % potenza di 2: passo finale = passo grosso / refine
camp(1).tFinal        = 600;
camp(1).modes         = {'nonlinear'};
camp(1).satCases      = true;

% Dominio ampio: delta rispetto al trim. alpha assoluto in [-80,80] deg,
% beta assoluto in [-145,155] deg (intorno ai limiti fisici/sensori).
camp(2).name          = 'GLOBALE';
camp(2).alphaRangeDeg = [-90 70];
camp(2).betaRangeDeg  = [-160 140];
camp(2).nAlpha        = 33;
camp(2).nBeta         = 61;
camp(2).refine        = 4;
camp(2).tFinal        = 600;
camp(2).modes         = {'nonlinear'};
camp(2).satCases      = [true false];

if opt.quick
    for ic = 1:numel(camp)
        camp(ic).nAlpha = 9;
        camp(ic).nBeta  = 9;
        camp(ic).refine = 2;
        camp(ic).tFinal = 60;
    end
    verify.factor  = 2;
    verify.nSample = 4;
    opt.saveMAT    = false;
end

%% ========================================================================
% 1. CONTROLLI SUL WORKSPACE E COSTRUZIONE DI cfg
% ========================================================================
reqVars = {'p','act','sensor','alpha_0','beta_0','u0_nom','x0','A_nom','B_nom','Cy'};
for k = 1:numel(reqVars)
    if ~exist(reqVars{k},'var')
        error('Variabile %s non presente nel workspace: eseguire prima il dataset.',reqVars{k});
    end
end
requiredP = {'J_alpha','l','eps_p','eps_y','J_y','J_z','m','I_b','c_alpha','c_beta','g'};
for k = 1:numel(requiredP)
    if ~isfield(p,requiredP{k})
        error('Parametro p.%s mancante: eseguire il dataset aggiornato.',requiredP{k});
    end
end

cfg.p        = p;
cfg.act      = act;
cfg.sensor   = sensor;
cfg.tol      = tol;
cfg.alpha_0  = alpha_0;
cfg.beta_0   = beta_0;
cfg.u0       = u0_nom(:);
cfg.x0       = x0(:);
cfg.A_lin    = double(A_nom);
cfg.B_lin    = double(B_nom);
cfg.Cy_lin   = double(Cy);
cfg.dt       = opt.dt;
cfg.verbose  = true;
cfg.recordTraj = false;
cfg.earlyStop  = opt.earlyStop;

if opt.f1NonNegative
    cfg.act.deltaF1_min = max(act.deltaF1_min,-u0_nom(1));
end

cfg.y0_sensor = [ ...
    -p.g*sin(alpha_0); ...
     cos(alpha_0)*cos(beta_0); ...
    -sin(beta_0)];

Jsensor = [ ...
    -p.g*cos(alpha_0),              0; ...
    -sin(alpha_0)*cos(beta_0), -cos(alpha_0)*sin(beta_0); ...
     0,                         -cos(beta_0)];
Vinv = diag(1./[sensor.acc.var; sensor.mag.var; sensor.mag.var]);
cfg.Hy = (Jsensor.'*Vinv*Jsensor) \ (Jsensor.'*Vinv);

if max(max(abs(Jsensor - cfg.Cy_lin(:,[1 3])))) > 1e-9
    warning('Jsensor diverso da Cy(:,[1 3]): controllare la definizione dei sensori.');
end

cfg.Am1 = [0 1; -act.wn1^2 -2*act.zeta1*act.wn1];
cfg.Bm1 = [0; act.wn1^2];
cfg.Cm1 = [1 0];
cfg.Am2 = [0 1; -act.wn2^2 -2*act.zeta2*act.wn2];
cfg.Bm2 = [0; act.wn2^2];
cfg.Cm2 = [1 0];

fprintf('\n============================================================\n');
fprintf('RAS V3 - trim: alpha0 = %.2f deg, beta0 = %.2f deg\n',rad2deg(alpha_0),rad2deg(beta_0));
fprintf('u0 = [%.4f; %.4f] N | dt = %.4f s | quick = %d\n',u0_nom(1),u0_nom(2),opt.dt,opt.quick);
fprintf('============================================================\n');

%% ========================================================================
% 2. CARICAMENTO CONTROLLORI, RIDUZIONE OPZIONALE, DUPLICATI
% ========================================================================
ctrls = loadctrls(opt.ctrlNames,opt.lqgRows);

if opt.muReduceOrder > 0
    for k = 1:numel(ctrls)
        if strcmp(ctrls(k).name,'mu-synthesis') && size(ctrls(k).A,1) > opt.muReduceOrder
            try
                Kfull = ss(ctrls(k).A,ctrls(k).B,ctrls(k).C,ctrls(k).D);
                Kred  = balred(Kfull,opt.muReduceOrder);
                Cr = ctrls(k);
                [Cr.A,Cr.B,Cr.C,Cr.D] = ssdata(Kred);
                [~,i1] = linearClosedLoop(ctrls(k),cfg);
                [~,i2] = linearClosedLoop(Cr,cfg);
                fprintf('\nRiduzione mu: ordine %d -> %d | max Re polo: %+.4f -> %+.4f\n', ...
                    size(ctrls(k).A,1),opt.muReduceOrder,i1.maxRe,i2.maxRe);
                if i2.maxRe < 0
                    ctrls(k) = Cr;
                    disp('  Riduzione accettata (loop lineare stabile). Verifica comunque i risultati.');
                else
                    disp('  Riduzione RIFIUTATA (loop lineare instabile): uso il controllore completo.');
                end
            catch ME
                fprintf('\nRiduzione mu fallita (%s): uso il controllore completo.\n',ME.message);
            end
        end
    end
end

dupOf = zeros(1,numel(ctrls));
for k = 2:numel(ctrls)
    for j = 1:k-1
        if dupOf(j) == 0 && sameController(ctrls(j),ctrls(k))
            dupOf(k) = j;
            break
        end
    end
end

fprintf('\n============================================================\n');
fprintf('CONTROLLORI CARICATI\n');
fprintf('============================================================\n');
for k = 1:numel(ctrls)
    fprintf('%-14s ordine = %3d, ingressi = %d, uscite = %2d, righe usate = %s', ...
        ctrls(k).name,size(ctrls(k).A,1),ctrls(k).nin,size(ctrls(k).C,1),mat2str(ctrls(k).rows));
    if dupOf(k) > 0
        fprintf('  [IDENTICO a %s]',ctrls(dupOf(k)).name);
    end
    fprintf('\n');
    if size(ctrls(k).C,1) > 2
        nm = ctrls(k).outNames;
        if iscell(nm) && ~all(cellfun(@isempty,nm))
            fprintf('    OutputName: %s\n',strjoin(nm(:).',', '));
        else
            fprintf('    ATTENZIONE: %d uscite senza OutputName. Verificare che le righe %s siano [F1;F2].\n', ...
                size(ctrls(k).C,1),mat2str(ctrls(k).rows));
        end
    end
end

%% ========================================================================
% 3. POLI DEL LOOP LINEARE (stabilita' locale e polo lento)
% ========================================================================
fprintf('\n============================================================\n');
fprintf('LOOP LINEARE (plant + attuatori + Pade'' 3o ordine + controllore)\n');
fprintf('============================================================\n');
linInfo = repmat(struct('maxRe',NaN,'slowRe',NaN,'tau',NaN),1,numel(ctrls));
for k = 1:numel(ctrls)
    [~,info] = linearClosedLoop(ctrls(k),cfg);
    linInfo(k) = info;
    fprintf('%-14s max Re(polo) = %+10.4f | polo lento = %+10.5f rad/s (tau = %8.1f s)\n', ...
        ctrls(k).name,info.maxRe,info.slowRe,info.tau);
    if info.maxRe >= 0
        fprintf('    ATTENZIONE: loop lineare NON asintoticamente stabile: la convergenza non e'' garantita neppure vicino al trim.\n');
    end
    if isfinite(info.tau) && info.tau > min([camp.tFinal])/12
        fprintf('    Nota: tau > Tf/12 -> il modo lento puo'' lasciare punti indeterminati (classe 2). Aumentare camp.tFinal.\n');
    end
end

%% ========================================================================
% 4. SELF-TEST: non lineare vs lineare su perturbazione piccola (0.2 deg)
% ========================================================================
if opt.selfTest
    fprintf('\n============================================================\n');
    fprintf('SELF-TEST lineare vs non lineare (d_alpha = 0.2 deg, d_beta = -0.2 deg, 10 s)\n');
    fprintf('============================================================\n');
    cfgT = cfg;
    cfgT.recordTraj = true;
    cfgT.verbose = false;
    cfgT.earlyStop = false;     % serve la traiettoria completa
    d0 = deg2rad(0.2);
    for k = 1:numel(ctrls)
        rn = simulateRAS(ctrls(k),'nonlinear',true,d0,-d0,10,cfgT,'selftest');
        rl = simulateRAS(ctrls(k),'linear',true,d0,-d0,10,cfgT,'selftest');
        Dn = rn.traj - cfg.x0;
        Dl = rl.traj;
        relErr = max(vecnorm(Dn-Dl,2,1))/max(vecnorm(Dl,2,1));
        fprintf('%-14s errore relativo max |NL - LIN| / |LIN| = %.4f', ctrls(k).name,relErr);
        if relErr > 0.05
            fprintf('   <-- ATTENZIONE: modello NL e linearizzazione non coerenti');
        end
        fprintf('\n');
    end
end

%% ========================================================================
% 5. CAMPAGNE DI SIMULAZIONE (griglia grossa + raffinamento del bordo)
% ========================================================================
runs = [];
for ic = 1:numel(camp)
    c = camp(ic);

    fprintf('\n############################################################\n');
    fprintf('CAMPAGNA %s: alpha [%g %g] deg (%d pt), beta [%g %g] deg (%d pt)\n', ...
        c.name,c.alphaRangeDeg(1),c.alphaRangeDeg(2),c.nAlpha, ...
        c.betaRangeDeg(1),c.betaRangeDeg(2),c.nBeta);
    fprintf('Raffinamento x%d -> reticolo fine %d x %d | passo finale %.3f x %.3f deg | Tf = %g s\n', ...
        c.refine,(c.nAlpha-1)*c.refine+1,(c.nBeta-1)*c.refine+1, ...
        diff(c.alphaRangeDeg)/((c.nAlpha-1)*c.refine), ...
        diff(c.betaRangeDeg)/((c.nBeta-1)*c.refine),c.tFinal);
    fprintf('############################################################\n');

    for im = 1:numel(c.modes)
        for is = 1:numel(c.satCases)
            useSat = logical(c.satCases(is));
            runIdx = zeros(1,numel(ctrls));
            for ik = 1:numel(ctrls)
                if opt.skipDuplicates && dupOf(ik) > 0 && runIdx(dupOf(ik)) > 0
                    res = runs(runIdx(dupOf(ik)));
                    res.name = ctrls(ik).name;
                    res.duplicateOf = ctrls(dupOf(ik)).name;
                    fprintf('\n[%s] %s identico a %s: risultato riutilizzato\n', ...
                        c.name,ctrls(ik).name,ctrls(dupOf(ik)).name);
                else
                    res = computeRASmap(ctrls(ik),c.modes{im},useSat,c,cfg,c.name);
                end
                if isempty(runs)
                    runs = res;
                else
                    runs(end+1) = res;
                end
                runIdx(ik) = numel(runs);
                if opt.saveMAT
                    save('RAS_v3_raw.mat','runs','camp','opt','tol');
                end
            end
        end
    end
end

%% ========================================================================
% 6. MAPPE
% ========================================================================
for ic = 1:numel(camp)
    for im = 1:numel(camp(ic).modes)
        for is = 1:numel(camp(ic).satCases)
            md = camp(ic).modes{im};
            sc = logical(camp(ic).satCases(is));
            sel = find(arrayfun(@(r) strcmp(r.campaign,camp(ic).name) && ...
                strcmpi(r.mode,md) && r.useSat == sc,runs));
            if isempty(sel)
                continue
            end
            if sc
                satTxt = 'con saturazione';
            else
                satTxt = 'senza saturazione';
            end
            ttl = sprintf('RAS %s - %s - %s',camp(ic).name,upper(md),satTxt);
            plotRunsFigure(runs(sel),ttl);
            plotOverlay(runs(sel),['Confronto bordi - ',ttl]);
        end
    end
end

%% ========================================================================
% 7. METRICHE
% ========================================================================
RASmetrics = buildMetrics(runs);
fprintf('\n============================================================\n');
fprintf('METRICHE RAS  (Inner = classe 1; Upper = classi 1+2)\n');
fprintf('============================================================\n');
disp(RASmetrics);
writetable(RASmetrics,'RAS_v3_metrics.csv');

for ir = 1:height(RASmetrics)
    if RASmetrics.TouchesDomainBoundary(ir) && RASmetrics.PctInner(ir) > 0
        fprintf('[%s | %s | %s | %s] la RAS tocca il bordo del dominio simulato: il confine vero e'' oltre.\n', ...
            RASmetrics.Campaign(ir),RASmetrics.Mode(ir),RASmetrics.Saturation(ir),RASmetrics.Controller(ir));
    end
    if RASmetrics.PctUndecided(ir) > 0
        fprintf('[%s | %s | %s | %s] %.2f %% dei punti e'' indeterminato (cl. 2): RAS vera tra %.2f %% e %.2f %%.\n', ...
            RASmetrics.Campaign(ir),RASmetrics.Mode(ir),RASmetrics.Saturation(ir),RASmetrics.Controller(ir), ...
            RASmetrics.PctUndecided(ir),RASmetrics.PctInner(ir),RASmetrics.PctUpper(ir));
    end
end

if opt.saveMAT
    save('RAS_v3_raw.mat','runs','camp','opt','tol','RASmetrics');
end

%% ========================================================================
% 8. VERIFICA FINALE
%    (a) punti indeterminati ricalcolati con orizzonte esteso;
%    (b) punti di bordo ricalcolati con dt/2 (la classe non deve dipendere
%        dal passo di integrazione).
% ========================================================================
verifyTable = table();
if verify.enable
    fprintf('\n============================================================\n');
    fprintf('VERIFICA FINALE (max %d punti per gruppo)\n',verify.nSample);
    fprintf('============================================================\n');
    rng(1);
    vrows = [];
    for ir = 1:numel(runs)
        r = runs(ir);
        if ~strcmpi(r.mode,'nonlinear') || ~r.useSat || ~isempty(r.duplicateOf)
            continue
        end
        k = find(strcmp({ctrls.name},r.name),1);
        [AA,BB] = meshgrid(r.aAxis,r.bAxis);
        dAall = deg2rad(AA(:)).';
        dBall = deg2rad(BB(:)).';
        M = r.clsMap;

        % (a) indeterminati -> orizzonte esteso
        idx = find(M(:) == 2).';
        if ~isempty(idx)
            idx = idx(randperm(numel(idx),min(verify.nSample,numel(idx))));
            cfgV = cfg;
            cfgV.verbose = false;
            rv = simulateRAS(ctrls(k),'nonlinear',true,dAall(idx),dBall(idx), ...
                verify.factor*r.tFinal,cfgV,[r.campaign,'-verifica']);
            fprintf('[%s] %-14s indeterminati, Tf x%g: %d/%d convergono, %d non convergono, %d ancora indeterminati\n', ...
                r.campaign,r.name,verify.factor,nnz(rv.cls == 1),numel(idx),nnz(rv.cls == 0),nnz(rv.cls == 2));
            vrow = struct('Campaign',string(r.campaign),'Controller',string(r.name), ...
                'Check',"indeterminati a Tf esteso",'nTested',numel(idx), ...
                'nConverged',nnz(rv.cls == 1),'nNotRAS',nnz(rv.cls == 0), ...
                'nUndecided',nnz(rv.cls == 2),'PctSameClass',NaN);
            if isempty(vrows)
                vrows = vrow;
            else
                vrows(end+1) = vrow;
            end
        end

        % (b) bordo con dt/2
        Bm = boundaryMask(M) & r.simMask;
        idx = find(Bm(:)).';
        if ~isempty(idx)
            idx = idx(randperm(numel(idx),min(verify.nSample,numel(idx))));
            cfgV = cfg;
            cfgV.verbose = false;
            cfgV.dt = cfg.dt/2;
            rv = simulateRAS(ctrls(k),'nonlinear',true,dAall(idx),dBall(idx), ...
                r.tFinal,cfgV,[r.campaign,'-verifica']);
            same = 100*mean(rv.cls == M(idx));
            fprintf('[%s] %-14s bordo con dt/2: %d punti, stessa classe nel %.1f %% dei casi\n', ...
                r.campaign,r.name,numel(idx),same);
            vrow = struct('Campaign',string(r.campaign),'Controller',string(r.name), ...
                'Check',"bordo con dt/2",'nTested',numel(idx), ...
                'nConverged',nnz(rv.cls == 1),'nNotRAS',nnz(rv.cls == 0), ...
                'nUndecided',nnz(rv.cls == 2),'PctSameClass',same);
            if isempty(vrows)
                vrows = vrow;
            else
                vrows(end+1) = vrow;
            end
        end
    end
    if ~isempty(vrows)
        verifyTable = struct2table(vrows);
        writetable(verifyTable,'RAS_v3_verify.csv');
    end
end

%% ========================================================================
% 9. SALVATAGGIO FINALE
% ========================================================================
if opt.saveMAT
    save('RAS_v3_results.mat','runs','camp','opt','tol','RASmetrics','verifyTable','linInfo');
    fprintf('\nSalvati: RAS_v3_raw.mat, RAS_v3_results.mat, RAS_v3_metrics.csv, RAS_v3_verify.csv\n');
end
disp('RAS V3 completata.');

%% ========================================================================
% FUNZIONI LOCALI
% ========================================================================
function res = computeRASmap(Ctl,mode,useSat,c,cfg,label)
    R = max(1,round(c.refine));
    if R ~= 2^round(log2(R))
        error('camp.refine deve essere una potenza di 2 (1, 2, 4, 8, ...).');
    end
    NA = (c.nAlpha-1)*R + 1;
    NB = (c.nBeta-1)*R + 1;
    aAx = linspace(c.alphaRangeDeg(1),c.alphaRangeDeg(2),NA);
    bAx = linspace(c.betaRangeDeg(1),c.betaRangeDeg(2),NB);

    S = nan(NB,NA);      % classi dei punti simulati (NaN = non simulato)
    F = nan(NB,NA);      % classi per riempimento delle celle uniformi
    acc.cmd1 = [];  acc.cmd2 = [];
    acc.sat1 = [];  acc.sat2 = [];
    acc.tConv = [];
    tStart = tic;

    % --- livello 0: griglia grossa ---
    [IB,IA] = ndgrid(1:R:NB,1:R:NA);
    [S,acc] = simulateLattice(S,IB(:),IA(:),aAx,bAx,Ctl,mode,useSat,c.tFinal,cfg, ...
        sprintf('%s grossa',label),acc);

    % --- livelli di raffinamento ---
    s = R;
    while true
        rowsC = 1:s:(NB-s);
        colsC = 1:s:(NA-s);
        c1 = S(rowsC,colsC);
        c2 = S(rowsC,colsC+s);
        c3 = S(rowsC+s,colsC);
        c4 = S(rowsC+s,colsC+s);
        known = ~(isnan(c1) | isnan(c2) | isnan(c3) | isnan(c4));
        uni   = known & (c1 == c2) & (c1 == c3) & (c1 == c4);
        mixed = known & ~uni;

        if s > 1
            [ir,jc] = find(uni);
            for q = 1:numel(ir)
                F(rowsC(ir(q))+(0:s),colsC(jc(q))+(0:s)) = c1(ir(q),jc(q));
            end
        end
        if s == 1 || ~any(mixed(:))
            break
        end

        h = s/2;
        need = false(NB,NA);
        [ir,jc] = find(mixed);
        for q = 1:numel(ir)
            need(rowsC(ir(q))+(0:h:s),colsC(jc(q))+(0:h:s)) = true;
        end
        need = need & isnan(S);
        [ib,ia] = find(need);
        if ~isempty(ib)
            [S,acc] = simulateLattice(S,ib,ia,aAx,bAx,Ctl,mode,useSat,c.tFinal,cfg, ...
                sprintf('%s raff. passo %d/%d',label,h,R),acc);
        end
        s = h;
    end

    % --- mappa finale: simulato dove disponibile, altrimenti riempimento ---
    clsMap = S;
    miss = isnan(clsMap);
    clsMap(miss) = F(miss);
    if any(isnan(clsMap(:)))
        warning('Punti del reticolo non assegnati: classificati come indeterminati (cl. 2).');
        clsMap(isnan(clsMap)) = 2;
    end

    res.campaign = label;
    res.name = Ctl.name;
    res.mode = mode;
    res.useSat = useSat;
    res.clsMap = clsMap;
    res.simMask = ~isnan(S);
    res.aAxis = aAx;
    res.bAxis = bAx;
    res.R = R;
    res.tFinal = c.tFinal;
    res.nSim = nnz(~isnan(S));
    res.maxCmd1 = acc.cmd1;
    res.maxCmd2 = acc.cmd2;
    res.satFrac1 = acc.sat1;
    res.satFrac2 = acc.sat2;
    res.tConv = acc.tConv;
    res.elapsed = toc(tStart);
    res.duplicateOf = '';

    nTot = NA*NB;
    fprintf('\n  >> [%s] %s | %s | sat=%d\n',label,upper(mode),Ctl.name,useSat);
    fprintf('     simulati %d punti su %d del reticolo fine (%.1f %%)\n', ...
        res.nSim,nTot,100*res.nSim/nTot);
    fprintf('     RAS (cl. 1) = %.2f %% | indeterminati (cl. 2) = %.2f %% | non convergono (cl. 0) = %.2f %%\n', ...
        100*nnz(clsMap == 1)/nTot,100*nnz(clsMap == 2)/nTot,100*nnz(clsMap == 0)/nTot);
    fprintf('     tempo totale %.1f s\n',res.elapsed);
end

function [S,acc] = simulateLattice(S,ib,ia,aAx,bAx,Ctl,mode,useSat,tFinal,cfg,lbl,acc)
% Simula i punti del reticolo (ib = indici beta/righe, ia = indici alpha/colonne).
    ib = ib(:);
    ia = ia(:);
    dA0 = deg2rad(aAx(ia));
    dB0 = deg2rad(bAx(ib));
    dA0 = dA0(:).';
    dB0 = dB0(:).';
    r = simulateRAS(Ctl,mode,useSat,dA0,dB0,tFinal,cfg,lbl);
    S(sub2ind(size(S),ib,ia)) = r.cls(:);
    acc.cmd1 = [acc.cmd1, r.maxCmd1];
    acc.cmd2 = [acc.cmd2, r.maxCmd2];
    acc.sat1 = [acc.sat1, r.satFrac1];
    acc.sat2 = [acc.sat2, r.satFrac2];
    acc.tConv = [acc.tConv, r.tConv];
end

function result = simulateRAS(Ctl,mode,useSat,dA0,dB0,tFinal,cfg,label)
% Simula N traiettorie in parallelo (vettorizzato) e le classifica:
%   1 = converge (stato completo a zero in modo permanente)
%   2 = indeterminata (decrescente ma non ancora convergente a tFinal)
%   0 = non converge (diverge, esce dai limiti, non decrescente)
    p   = cfg.p;
    act = cfg.act;
    sensor = cfg.sensor;
    tol = cfg.tol;
    dt  = cfg.dt;

    isNL  = strcmpi(mode,'nonlinear');
    N     = numel(dA0);
    nIter = round(tFinal/dt);
    tailStart = tFinal*(1-tol.tailFraction);
    tailMid   = tailStart + 0.5*(tFinal-tailStart);
    maxAlpha  = deg2rad(tol.maxAlphaDeg);
    maxBeta   = deg2rad(tol.maxBetaDeg);

    earlyStop = cfg.earlyStop;
    aConv    = deg2rad(tol.convAngleDeg);
    rConv    = deg2rad(tol.convRateDeg);
    chkEvery = max(1,round(tol.checkEvery/dt));
    pkEvery  = max(1,round(0.02/dt));

    if cfg.verbose
        fprintf('\n------------------------------------------------------------\n');
        fprintf('RAS [%s] %s | %s | sat=%d | N=%d | Tf=%g s | dt=%g s\n', ...
            label,upper(mode),Ctl.name,useSat,N,tFinal,dt);
    end
    tic;

    % --- controllore: propagazione esatta a ingresso costante nel passo ---
    nk = size(Ctl.A,1);
    Cc = Ctl.C(Ctl.rows,:);
    Dc = Ctl.D(Ctl.rows,:);
    [Phi,Gam]   = exactZOH(Ctl.A,Ctl.B,dt);
    [PhiH,GamH] = exactZOH(Ctl.A,Ctl.B,dt/2);

    % --- plant lineare: propagazione esatta ---
    if ~isNL
        [PhiLin,GamLin] = exactZOH(cfg.A_lin,cfg.B_lin,dt);
    end

    Am1 = cfg.Am1; Bm1 = cfg.Bm1; Cm1 = cfg.Cm1;
    Am2 = cfg.Am2; Bm2 = cfg.Bm2; Cm2 = cfg.Cm2;

    % --- delay come buffer circolare ---
    dS1 = round(act.td1/dt);
    dS2 = round(act.td2/dt);
    if abs(dS1*dt-act.td1) > 1e-9 || abs(dS2*dt-act.td2) > 1e-9
        warning('dt non divide esattamente il transport delay: il delay viene quantizzato.');
    end
    dS1 = max(dS1,1);
    dS2 = max(dS2,1);
    dBuf1 = zeros(dS1,N);
    dBuf2 = zeros(dS2,N);
    ptr1 = 1;
    ptr2 = 1;

    accEvery = max(1,round(sensor.acc.Ts/dt));
    magEvery = max(1,round(sensor.mag.Ts/dt));

    % --- stati ---
    if isNL
        X = repmat(cfg.x0,1,N);
        X(1,:) = cfg.alpha_0 + dA0;
        X(3,:) = cfg.beta_0  + dB0;
        xeq = cfg.x0;
    else
        X = zeros(4,N);
        X(1,:) = dA0;
        X(3,:) = dB0;
        xeq = zeros(4,1);
    end
    Xk  = zeros(nk,N);
    Xm1 = zeros(2,N);
    Xm2 = zeros(2,N);
    yHold = zeros(3,N);     % misure delta-y con ZOH separati acc / mag

    alive     = true(1,N);
    discarded = false(1,N);
    converged = false(1,N);
    tConv     = nan(1,N);
    calmTime  = zeros(1,N);
    peakXk    = zeros(1,N);
    peakXm    = zeros(1,N);
    peakBuf   = zeros(1,N);
    tailNorm1 = zeros(1,N);
    tailNorm2 = zeros(1,N);
    maxCmd1 = zeros(1,N);
    maxCmd2 = zeros(1,N);
    satCnt1 = zeros(1,N);
    satCnt2 = zeros(1,N);
    stepsRun = zeros(1,N);

    trajLog = [];
    if cfg.recordTraj
        trajLog = zeros(4,nIter);
    end

    for istep = 1:nIter
        traj = find(alive);
        if isempty(traj)
            break
        end
        n  = numel(traj);
        Xt = X(:,traj);
        stepsRun(traj) = stepsRun(traj) + 1;

        % -------------------------------------------------------------
        % MISURE: ogni sensore ha il proprio ZOH
        % -------------------------------------------------------------
        updAcc = (mod(istep-1,accEvery) == 0);
        updMag = (mod(istep-1,magEvery) == 0);
        if updAcc || updMag
            if isNL
                ynow = [ ...
                    -p.g*sin(Xt(1,:)); ...
                     cos(Xt(1,:)).*cos(Xt(3,:)); ...
                    -sin(Xt(3,:))] - cfg.y0_sensor;
            else
                ynow = cfg.Cy_lin*Xt;
            end
            if updAcc
                yHold(1,traj) = ynow(1,:);
            end
            if updMag
                yHold(2:3,traj) = ynow(2:3,:);
            end
        end
        deltaY = yHold(:,traj);

        % -------------------------------------------------------------
        % INGRESSO AL CONTROLLORE
        % -------------------------------------------------------------
        switch Ctl.nin
            case 5
                uk = [zeros(2,n); deltaY];          % [r_alpha; r_beta; delta_y]
            case 3
                uk = deltaY;                        % LQG 1-DOF
            case 2
                uk = -(cfg.Hy*deltaY);              % errore angolare stimato (WLS)
            otherwise
                error('Numero di ingressi non riconosciuto per %s.',Ctl.name);
        end

        % -------------------------------------------------------------
        % PLANT (forze ritardate)
        % -------------------------------------------------------------
        dF1 = dBuf1(ptr1,traj);
        dF2 = dBuf2(ptr2,traj);
        if isNL
            Fp = [cfg.u0(1)+dF1; cfg.u0(2)+dF2];
            X(:,traj) = rk4Step(Xt,Fp,p,dt);
        else
            Fp = [dF1; dF2];
            X(:,traj) = PhiLin*Xt + GamLin*Fp;
        end
        if cfg.recordTraj
            trajLog(:,istep) = X(:,1);
        end

        % -------------------------------------------------------------
        % CONTROLLORE + ATTUATORI (RK4 a 3 punti)
        % -------------------------------------------------------------
        Xk0 = Xk(:,traj);
        if nk > 0
            XkH = PhiH*Xk0 + GamH*uk;
            XkE = Phi*Xk0  + Gam*uk;
        else
            XkH = Xk0;
            XkE = Xk0;
        end
        [c1a,c2a,r1a,r2a] = ctrlCommand(Xk0,uk,Cc,Dc,act,useSat);
        [c1b,c2b]         = ctrlCommand(XkH,uk,Cc,Dc,act,useSat);
        [c1c,c2c]         = ctrlCommand(XkE,uk,Cc,Dc,act,useSat);

        xm = Xm1(:,traj);
        k1 = Am1*xm + Bm1*c1a;
        k2 = Am1*(xm+0.5*dt*k1) + Bm1*c1b;
        k3 = Am1*(xm+0.5*dt*k2) + Bm1*c1b;
        k4 = Am1*(xm+dt*k3) + Bm1*c1c;
        Xm1(:,traj) = xm + (dt/6)*(k1+2*k2+2*k3+k4);

        xm = Xm2(:,traj);
        k1 = Am2*xm + Bm2*c2a;
        k2 = Am2*(xm+0.5*dt*k1) + Bm2*c2b;
        k3 = Am2*(xm+0.5*dt*k2) + Bm2*c2b;
        k4 = Am2*(xm+dt*k3) + Bm2*c2c;
        Xm2(:,traj) = xm + (dt/6)*(k1+2*k2+2*k3+k4);

        if nk > 0
            Xk(:,traj) = XkE;
        end

        maxCmd1(traj) = max(maxCmd1(traj),max([abs(c1a);abs(c1b);abs(c1c)],[],1));
        maxCmd2(traj) = max(maxCmd2(traj),max([abs(c2a);abs(c2b);abs(c2c)],[],1));
        satCnt1(traj) = satCnt1(traj) + double(r1a < act.deltaF1_min | r1a > act.deltaF1_max);
        satCnt2(traj) = satCnt2(traj) + double(r2a < act.deltaF2_min | r2a > act.deltaF2_max);

        % -------------------------------------------------------------
        % DELAY BUFFER
        % -------------------------------------------------------------
        dBuf1(ptr1,traj) = Cm1*Xm1(:,traj);
        dBuf2(ptr2,traj) = Cm2*Xm2(:,traj);
        ptr1 = ptr1 + 1;
        if ptr1 > dS1
            ptr1 = 1;
        end
        ptr2 = ptr2 + 1;
        if ptr2 > dS2
            ptr2 = 1;
        end

        % -------------------------------------------------------------
        % SCARTO (non finito o fuori dai limiti fisici)
        % -------------------------------------------------------------
        Xn = X(:,traj);
        if isNL
            aPh = Xn(1,:);
            bPh = Xn(3,:);
        else
            aPh = cfg.alpha_0 + Xn(1,:);
            bPh = cfg.beta_0  + Xn(3,:);
        end

        fin = all(isfinite(Xn),1) & all(isfinite(Xm1(:,traj)),1) & all(isfinite(Xm2(:,traj)),1);
        if nk > 0
            fin = fin & all(isfinite(Xk(:,traj)),1);
        end
        bad = ~fin | abs(aPh) >= maxAlpha | abs(bPh) >= maxBeta;
        if any(bad)
            discarded(traj(bad)) = true;
            alive(traj(bad)) = false;
        end

        % -------------------------------------------------------------
        % PICCHI DEI BLOCCHI NON-IMPIANTO (scala per il criterio relativo)
        % -------------------------------------------------------------
        if earlyStop && mod(istep,pkEvery) == 0
            tr = find(alive);
            if ~isempty(tr)
                if nk > 0
                    peakXk(tr) = max(peakXk(tr),max(abs(Xk(:,tr)),[],1));
                end
                peakXm(tr)  = max(peakXm(tr),max(abs([Xm1(:,tr);Xm2(:,tr)]),[],1));
                peakBuf(tr) = max(peakBuf(tr),max(abs([dBuf1(:,tr);dBuf2(:,tr)]),[],1));
            end
        end

        % -------------------------------------------------------------
        % CONVERGENZA SULLO STATO COMPLETO
        % -------------------------------------------------------------
        if earlyStop && mod(istep,chkEvery) == 0
            tr = find(alive);
            if ~isempty(tr)
                dX  = X(:,tr) - xeq;
                okP = all(abs(dX([1 3],:)) < aConv,1) & all(abs(dX([2 4],:)) < rConv,1);
                nXm = max(abs([Xm1(:,tr);Xm2(:,tr)]),[],1);
                nBf = max(abs([dBuf1(:,tr);dBuf2(:,tr)]),[],1);
                peakXm(tr)  = max(peakXm(tr),nXm);
                peakBuf(tr) = max(peakBuf(tr),nBf);
                okB = (nXm <= tol.relBlock*peakXm(tr)) & (nBf <= tol.relBlock*peakBuf(tr));
                if nk > 0
                    nXk = max(abs(Xk(:,tr)),[],1);
                    peakXk(tr) = max(peakXk(tr),nXk);
                    okB = okB & (nXk <= tol.relBlock*peakXk(tr));
                end
                calm = okP & okB;
                calmTime(tr) = (calmTime(tr) + tol.checkEvery).*calm;
                done = calmTime(tr) >= tol.calmWindow;
                if any(done)
                    converged(tr(done)) = true;
                    tConv(tr(done)) = istep*dt;
                    alive(tr(done)) = false;
                end
            end
        end

        % -------------------------------------------------------------
        % STATISTICHE DI CODA 
        % -------------------------------------------------------------
        tNext = istep*dt;
        if tNext >= tailStart
            tr = find(alive);
            if ~isempty(tr)
                dX  = X(:,tr) - xeq;
                nrm = vecnorm(dX,2,1);
                if tNext <= tailMid
                    tailNorm1(tr) = max(tailNorm1(tr),nrm);
                else
                    tailNorm2(tr) = max(tailNorm2(tr),nrm);
                end
            end
        end

        if cfg.verbose && mod(istep,max(1,round(nIter/10))) == 0
            fprintf('  %5.1f %% | in corso: %d / %d\n',100*istep/nIter,nnz(alive),N);
        end
    end

    % --- classificazione finale ---
    valid = ~discarded;
    cls = zeros(1,N);
    undecided = valid & ~converged & (tailNorm1 > 0) & (tailNorm2 < tol.decayFactor*tailNorm1);
    cls(undecided) = 2;
    cls(converged) = 1;

    elapsed = toc;
    if cfg.verbose
        fprintf('  classi: 0 non converge=%d | 1 converge (RAS)=%d | 2 indeterminata=%d\n', ...
            nnz(cls == 0),nnz(cls == 1),nnz(cls == 2));
        fprintf('  RAS (cl. 1): %d / %d (%.1f %%) | massimo (cl. 1+2): %.1f %%\n', ...
            nnz(cls == 1),N,100*nnz(cls == 1)/N,100*nnz(cls >= 1)/N);
        fprintf('  comando max: F1 = %.3f N, F2 = %.3f N | tempo simulato in saturazione: %.1f %% / %.1f %%\n', ...
            max(maxCmd1),max(maxCmd2), ...
            100*mean(satCnt1./max(stepsRun,1)),100*mean(satCnt2./max(stepsRun,1)));
        fprintf('  completato in %.1f s\n',elapsed);
    end

    result.cls = cls;
    result.tConv = tConv;
    result.maxCmd1 = maxCmd1;
    result.maxCmd2 = maxCmd2;
    result.satFrac1 = satCnt1./max(stepsRun,1);
    result.satFrac2 = satCnt2./max(stepsRun,1);
    result.elapsed = elapsed;
    result.traj = trajLog;
end

function B = boundaryMask(M)
% Punti del reticolo con almeno un vicino (4-connesso) di classe diversa.
    B = false(size(M));
    d1 = M(1:end-1,:) ~= M(2:end,:);
    d2 = M(:,1:end-1) ~= M(:,2:end);
    B(1:end-1,:) = B(1:end-1,:) | d1;
    B(2:end,:)   = B(2:end,:)   | d1;
    B(:,1:end-1) = B(:,1:end-1) | d2;
    B(:,2:end)   = B(:,2:end)   | d2;
end

function [CL,info] = linearClosedLoop(Ctl,cfg)
% Loop lineare: plant (delta) + attuatori 2o ordine + Pade' 3o ordine + K.
    act = cfg.act;
    P = ss(cfg.A_lin,cfg.B_lin,cfg.Cy_lin,zeros(3,2));
    [n1,d1] = pade(act.td1,3);
    [n2,d2] = pade(act.td2,3);
    Act1 = ss(tf(n1,d1))*ss(cfg.Am1,cfg.Bm1,cfg.Cm1,0);
    Act2 = ss(tf(n2,d2))*ss(cfg.Am2,cfg.Bm2,cfg.Cm2,0);
    G = P*blkdiag(Act1,Act2);

    Cc = Ctl.C(Ctl.rows,:);
    Dc = Ctl.D(Ctl.rows,:);
    switch Ctl.nin
        case 5
            Kc = ss(Ctl.A,Ctl.B(:,3:5),Cc,Dc(:,3:5));
            CL = feedback(G,Kc,+1);
        case 3
            Kc = ss(Ctl.A,Ctl.B,Cc,Dc);
            CL = feedback(G,Kc,+1);
        case 2
            Kc = ss(Ctl.A,Ctl.B,Cc,Dc)*cfg.Hy;
            CL = feedback(G,Kc,-1);
        otherwise
            error('Numero di ingressi non riconosciuto per %s.',Ctl.name);
    end

    pl = pole(CL);
    info.maxRe = max(real(pl));
    stab = pl(real(pl) < 0);
    if isempty(stab)
        info.slowRe = NaN;
        info.tau = NaN;
    else
        [~,i] = min(abs(real(stab)));
        info.slowRe = real(stab(i));
        info.tau = 1/abs(info.slowRe);
    end
end

function same = sameController(C1,C2)
% true se i due controllori sono identici (stessa realizzazione o stessa FdT).
    same = false;
    if C1.nin ~= C2.nin || any(size(C1.A) ~= size(C2.A)) || any(size(C1.C) ~= size(C2.C))
        return
    end
    mx = @(M) max([0; abs(M(:))]);
    sc = 1 + max([mx(C1.A),mx(C1.B),mx(C1.C),mx(C1.D)]);
    d  = max([mx(C1.A-C2.A),mx(C1.B-C2.B),mx(C1.C-C2.C),mx(C1.D-C2.D)]);
    if d < 1e-9*sc
        same = true;
        return
    end
    try
        K1 = ss(C1.A,C1.B,C1.C(C1.rows,:),C1.D(C1.rows,:));
        K2 = ss(C2.A,C2.B,C2.C(C2.rows,:),C2.D(C2.rows,:));
        n1 = norm(K1,inf);
        if isfinite(n1) && n1 > 0
            same = norm(K1-K2,inf)/n1 < 1e-6;
        end
    catch
        same = false;
    end
end

function plotRunsFigure(runsSel,figName)
    % classi: 0 rosso (non converge), 1 verde (RAS), 2 giallo (indeterminata)
    cmap = [0.80 0.15 0.15; 0.10 0.55 0.25; 0.98 0.85 0.20];
    figure('Name',figName,'Color','w');
    tl = tiledlayout('flow','TileSpacing','compact');
    title(tl,figName,'Interpreter','none');
    for k = 1:numel(runsSel)
        r = runsSel(k);
        nexttile;
        ax = gca;
        imagesc(r.aAxis,r.bAxis,r.clsMap,[0 2]);
        axis xy;
        hold on; grid on; box on;
        colormap(ax,cmap);
        A = double(r.clsMap == 1);
        if any(A(:)) && ~all(A(:))
            contour(r.aAxis,r.bAxis,A,[0.5 0.5],'k','LineWidth',1.5);
        end
        plot(0,0,'k+','MarkerSize',10,'LineWidth',2);
        xlabel('$\Delta\alpha(0)\;[^{\circ}]$','Interpreter','latex');
        ylabel('$\Delta\beta(0)\;[^{\circ}]$','Interpreter','latex');
        t = r.name;
        if ~isempty(r.duplicateOf)
            t = [t,' (= ',r.duplicateOf,')'];
        end
        title(sprintf('%s - RAS %.1f %% (indet. %.1f %%)',t, ...
            100*mean(r.clsMap(:) == 1),100*mean(r.clsMap(:) == 2)),'Interpreter','none');
        xlim(r.aAxis([1 end]));
        ylim(r.bAxis([1 end]));
    end
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Ticks = 0:2;
    cb.TickLabels = {'0 non converge','1 converge (RAS)','2 indeterminata'};
end

function plotOverlay(runsSel,figName)
    figure('Name',figName,'Color','w');
    hold on; grid on; box on;
    cols = lines(numel(runsSel));
    for k = 1:numel(runsSel)
        r = runsSel(k);
        if ~isempty(r.duplicateOf)
            continue
        end
        A = double(r.clsMap == 1);
        if any(A(:)) && ~all(A(:))
            [~,h] = contour(r.aAxis,r.bAxis,A,[0.5 0.5],'LineWidth',2,'DisplayName',r.name);
            h.LineColor = cols(k,:);
        elseif all(A(:))
            rectangle('Position',[r.aAxis(1) r.bAxis(1) r.aAxis(end)-r.aAxis(1) r.bAxis(end)-r.bAxis(1)], ...
                'EdgeColor',cols(k,:),'LineStyle','--','LineWidth',1.5);
            plot(NaN,NaN,'--','Color',cols(k,:),'LineWidth',1.5, ...
                'DisplayName',[r.name,' (tutto il dominio)']);
        end
    end
    r1 = runsSel(1);
    plot(0,0,'k+','MarkerSize',10,'LineWidth',2,'HandleVisibility','off');
    xlabel('$\Delta\alpha(0)\;[^{\circ}]$','Interpreter','latex');
    ylabel('$\Delta\beta(0)\;[^{\circ}]$','Interpreter','latex');
    title(figName,'Interpreter','none');
    xlim(r1.aAxis([1 end]));
    ylim(r1.bAxis([1 end]));
    legend('Location','bestoutside','Interpreter','none');
end

function T = buildMetrics(runs)
    rows = [];
    for ir = 1:numel(runs)
        r = runs(ir);
        M = r.clsMap;
        N = numel(M);
        dAs = mean(diff(r.aAxis));
        dBs = mean(diff(r.bAxis));
        cellA = dAs*dBs;
        inner = (M == 1);
        upper = (M >= 1);

        row = struct();
        row.Campaign = string(r.campaign);
        row.Mode = string(upper(r.mode));
        if r.useSat
            row.Saturation = "sat";
        else
            row.Saturation = "nosat";
        end
        row.Controller = string(r.name);
        row.DuplicateOf = string(r.duplicateOf);
        row.N = N;
        row.NSimulated = r.nSim;
        row.StepAlphaDeg = dAs;
        row.StepBetaDeg = dBs;
        row.PctInner = 100*nnz(inner)/N;
        row.PctUndecided = 100*nnz(M == 2)/N;
        row.PctUpper = 100*nnz(upper)/N;
        row.PctOutside = 100*nnz(M == 0)/N;
        row.AreaInnerDeg2 = nnz(inner)*cellA;
        row.AreaUpperDeg2 = nnz(upper)*cellA;
        row.DeltaAlphaMinDeg = NaN;
        row.DeltaAlphaMaxDeg = NaN;
        row.DeltaBetaMinDeg = NaN;
        row.DeltaBetaMaxDeg = NaN;
        row.TouchesDomainBoundary = false;
        if any(inner(:))
            [ib,ia] = find(inner);
            row.DeltaAlphaMinDeg = min(r.aAxis(ia));
            row.DeltaAlphaMaxDeg = max(r.aAxis(ia));
            row.DeltaBetaMinDeg = min(r.bAxis(ib));
            row.DeltaBetaMaxDeg = max(r.bAxis(ib));
            edgeMask = [inner(1,:),inner(end,:),inner(:,1).',inner(:,end).'];
            row.TouchesDomainBoundary = any(edgeMask);
        end
        row.MaxCmdF1 = max(r.maxCmd1);
        row.MaxCmdF2 = max(r.maxCmd2);
        row.MeanSatFracF1 = mean(r.satFrac1);
        row.MeanSatFracF2 = mean(r.satFrac2);
        row.MedianTConv_s = median(r.tConv,'omitnan');

        if isempty(rows)
            rows = row;
        else
            rows(end+1) = row;
        end
    end
    if isempty(rows)
        T = table();
    else
        T = struct2table(rows);
    end
end

function ctrls = loadctrls(reqNames,lqgRows)
    ctrls = struct('name',{},'family',{},'A',{},'B',{},'C',{},'D',{}, ...
                   'nin',{},'rows',{},'outNames',{});

    if isfile('LQG_Controllers.mat')
        SL = load('LQG_Controllers.mat');
    else
        SL = struct;
    end
    if ismember('LQG1',reqNames)
        K = fetchSystem('K_lqg_reg',SL);
        ctrls(end+1) = makectrl('LQG 1-DOF','LQG',K,lqgRows);
    end
    if ismember('LQG2',reqNames)
        K = fetchSystem('K_lqg_2dof',SL);
        ctrls(end+1) = makectrl('LQG 2-DOF','LQG',K,lqgRows);
    end
    if ismember('LQGI',reqNames)
        K = fetchSystem('K_lqg_int',SL);
        ctrls(end+1) = makectrl('LQGI','LQG',K,lqgRows);
    end

    if isfile('HINF_controllers.mat')
        SH = load('HINF_controllers.mat');
    else
        SH = struct;
    end
    if ismember('mixsyn',reqNames)
        K = fetchSystem('K_mix',SH);
        ctrls(end+1) = makectrl('mixsyn','HINF',K,lqgRows);
    end
    if ismember('hinfsyn',reqNames)
        K = fetchSystem('K_hinfsyn',SH);
        ctrls(end+1) = makectrl('hinfsyn','HINF',K,lqgRows);
    end
    if ismember('PID+comp',reqNames)
        K = fetchSystem('K_pidcomp',SH);
        ctrls(end+1) = makectrl('PID+comp','HINF',K,lqgRows);
    end

    if ismember('mu-synthesis',reqNames)
        if isfile('MU_controller.mat')
            SM = load('MU_controller.mat');
        else
            SM = struct;
        end
        K = fetchSystem('K_mu',SM);
        ctrls(end+1) = makectrl('mu-synthesis','HINF',K,lqgRows);
    end

    if isempty(ctrls)
        error(['Nessun controllore richiesto trovato. ', ...
            'Controllare opt.ctrlNames e i file .mat.']);
    end
end

function K = fetchSystem(varName,S)
    if isfield(S,varName)
        K = S.(varName);
        return
    end
    if evalin('base',sprintf('exist(''%s'',''var'')',varName))
        K = evalin('base',varName);
        return
    end
    error('Controllore %s non trovato ne'' nel file .mat ne'' nel workspace.',varName);
end

function Ctl = makectrl(name,family,K,lqgRows)
    K_ss = ss(K);
    [A,B,C,D] = ssdata(K_ss);
    Ctl.name = name;
    Ctl.family = family;
    Ctl.A = double(A);
    Ctl.B = double(B);
    Ctl.C = double(C);
    Ctl.D = double(D);
    Ctl.nin = size(Ctl.B,2);
    if size(Ctl.C,1) == 2
        Ctl.rows = [1 2];
    else
        Ctl.rows = lqgRows;
    end
    Ctl.outNames = K_ss.OutputName;
end

function [Phi,Gamma] = exactZOH(A,B,h)
    n = size(A,1);
    m = size(B,2);
    if n == 0
        Phi = zeros(0);
        Gamma = zeros(0,m);
        return
    end
    E = expm([A,B;zeros(m,n+m)]*h);
    Phi = E(1:n,1:n);
    Gamma = E(1:n,n+(1:m));
end

function [cmd1,cmd2,raw1,raw2] = ctrlCommand(Xk,uk,Cc,Dc,act,useSat)
    raw = Dc*uk + Cc*Xk;
    raw1 = raw(1,:);
    raw2 = raw(2,:);
    if useSat
        cmd1 = min(max(raw1,act.deltaF1_min),act.deltaF1_max);
        cmd2 = min(max(raw2,act.deltaF2_min),act.deltaF2_max);
    else
        cmd1 = raw1;
        cmd2 = raw2;
    end
end

function dX = helicopterDerivativeVectorized(X,F,p)
    alpha = X(1,:);
    alphaDot = X(2,:);
    beta = X(3,:);
    betaDot = X(4,:);
    F1 = F(1,:);
    F2 = F(2,:);

    Jbeta = p.J_y.*sin(alpha).^2 + ...
        (p.J_z+p.m*p.l^2).*cos(alpha).^2 + p.I_b;

    tauAlpha = p.l.*(cos(beta).*F1 + p.eps_p.*sin(beta).*F2) ...
        - p.c_alpha.*alphaDot - p.m*p.g*p.l.*sin(alpha);

    tauBeta = p.l.*(p.eps_y.*sin(alpha).*F1 + cos(alpha).*F2) ...
        - p.c_beta.*betaDot;

    dX = [alphaDot; tauAlpha./p.J_alpha; betaDot; tauBeta./Jbeta];
end

function Xnext = rk4Step(X,F,p,h)
    k1 = helicopterDerivativeVectorized(X,F,p);
    k2 = helicopterDerivativeVectorized(X+0.5*h*k1,F,p);
    k3 = helicopterDerivativeVectorized(X+0.5*h*k2,F,p);
    k4 = helicopterDerivativeVectorized(X+h*k3,F,p);
    Xnext = X + h*(k1+2*k2+2*k3+k4)/6;
end
