%% =========================================================================
%  CONFRONTO_RISULTATI - famiglia LQG
%  LQG 2DOF | H2-LQG | H2-soft | LQG con integratore
% =========================================================================
clc; close all;

%% 0. CONFIGURAZIONE
cfg.tRef     = ref.alpha.time;
cfg.tDist    = aero.alpha.time;
cfg.refAlpha = ref.alpha.final;            % [rad]
cfg.refBeta  = ref.beta.final;             % [rad]
cfg.bandTrk  = 0.02;                       % banda di assestamento (2%)
cfg.bandDist = deg2rad(0.05);              % banda di rientro dopo il disturbo [rad]
cfg.tail     = 5;                          % finestra finale per l'offset [s]
cfg.dt       = 0.005;
cfg.limits   = [act.deltaF1_min act.deltaF1_max; ...
                act.deltaF2_min act.deltaF2_max];
cfg.dAmp     = [aero.alpha.amplitude; aero.beta.amplitude];   % ampiezze disturbi [N*m]
cfg.outDir   = ['figure_confronto_' datestr(now,'yyyymmdd_HHMMSS')];
mkdir(cfg.outDir);
diary(fullfile(cfg.outDir, 'log_console.txt'));   % salva anche tutto cio' che viene stampato
diary on

%% CASI: etichetta, nome dei To Workspace, versione
defs = struct('label',{'LQG 2DOF','H2-LQG','H2-soft','LQG integratore'}, ...
              'tag',  {'LQG2',    'LQG_H2','LQGsoft','LQG_int'});
vars = struct('suffix',{'','_nl'}, 'group',{'Lineare','Non lineare'}, 'lab',{'',' (NL)'});

cases = struct('label',{},'key',{},'group',{},'nameA',{},'nameB',{},'nameU',{},'isNL',{});
fprintf('Ricerca dei segnali nel workspace...\n');
for v = 1:numel(vars)
    for d = 1:numel(defs)
        nA = ['alpha_' defs(d).tag vars(v).suffix];
        nB = ['beta_'  defs(d).tag vars(v).suffix];
        nU = ['u_'     defs(d).tag vars(v).suffix];
        if hasVar(nA) && hasVar(nB)
            cases(end+1) = struct('label',[defs(d).label vars(v).lab], 'key',defs(d).label, ...
                                  'group',vars(v).group, 'nameA',nA, 'nameB',nB, ...
                                  'nameU',nU, 'isNL',v==2);
            fprintf('  trovato: %s, %s\n', nA, nB);
        else
            fprintf('  (salto) %s / %s non trovati\n', nA, nB);
        end
    end
end
assert(~isempty(cases), 'Nessun segnale trovato: lancia prima il Simulink.');

% coppie da confrontare punto a punto (equivalenza nel tempo)
pairs = { 'LQG 2DOF','H2-LQG'; 'LQG 2DOF (NL)','H2-LQG (NL)' };

%% ========================================================================
%  DOMINIO DEL TEMPO
%% ========================================================================
nC  = numel(cases);
raw = cell(nC,1);
tEnd = inf;
for k = 1:nC
    [ta,a] = getsig(cases(k).nameA);
    [tb,b] = getsig(cases(k).nameB);
    r = struct('ta',ta,'a',a(:,1),'tb',tb,'b',b(:,1),'tu',[],'u',[]);
    if hasVar(cases(k).nameU)
        [tu,u] = getsig(cases(k).nameU);
        r.tu = tu;  r.u = u(:,1:2);
    end
    raw{k} = r;
    tEnd = min([tEnd, ta(end), tb(end)]);
end
cfg.tStop = tEnd;
tc = (0:cfg.dt:tEnd).';

R = repmat(struct('a',[],'b',[],'u',[]), nC, 1);
for k = 1:nC
    a = raw{k}.a;  b = raw{k}.b;
    % Segnali assoluti (tipico del non lineare): si sottrae il trim
    trim = [0 0];
    if exist('alpha_0','var') && exist('beta_0','var') && ...
       abs(a(1) - alpha_0) < 0.25*abs(alpha_0) && abs(b(1) - beta_0) < 0.25*abs(beta_0)
        trim = [alpha_0 beta_0];
        fprintf('  %s: segnali assoluti, sottraggo il trim [%.3f %.3f] rad\n', cases(k).label, trim);
    end
    R(k).a = interp1(raw{k}.ta, a - trim(1), tc, 'linear');
    R(k).b = interp1(raw{k}.tb, b - trim(2), tc, 'linear');
    if ~isempty(raw{k}.u)
        R(k).u = interp1(raw{k}.tu, raw{k}.u, tc, 'linear');
    end
end

rows = cell(nC,1);
for k = 1:nC
    rows{k} = metrics(tc, R(k).a, R(k).b, R(k).u, cfg);
    rows{k}.Caso = string(cases(k).label);
end
Tt = struct2table([rows{:}]);
Tt = movevars(Tt, 'Caso', 'Before', 1);
disp('=== METRICHE NEL TEMPO ===');
disp(Tt);

% --- equivalenza punto a punto
fprintf('\n=== EQUIVALENZA PUNTO A PUNTO ===\n');
Te = table();
for j= 1:size(pairs,1)
    i1 = find(strcmp({cases.label}, pairs{j,1}), 1);
    i2 = find(strcmp({cases.label}, pairs{j,2}), 1);
    if isempty(i1) || isempty(i2), continue; end
    dA = max(abs(R(i1).a-R(i2).a));   dB = max(abs(R(i1).b-R(i2).b));
    fprintf('%s vs %s: max|d alpha| = %.3e rad, max|d beta| = %.3e rad\n', ...
        pairs{j,1}, pairs{j,2}, dA, dB);
    Te = [Te; table(string(pairs{j,1}), string(pairs{j,2}), dA, dB, ...
            'VariableNames',{'Caso_1','Caso_2','MaxDiff_alpha_rad','MaxDiff_beta_rad'})]; 
end

% --- scarto lineare / non lineare per ogni controllore
Tn = table();
for d = 1:numel(defs)
    iL = find(strcmp({cases.label}, defs(d).label), 1);
    iN = find(strcmp({cases.label}, [defs(d).label ' (NL)']), 1);
    if isempty(iL) || isempty(iN), continue; end
    row = table(string(defs(d).label), ...
        rad2deg(max(abs(R(iL).a - R(iN).a))), rad2deg(max(abs(R(iL).b - R(iN).b))), ...
        rad2deg(rms(R(iL).a - R(iN).a)),      rad2deg(rms(R(iL).b - R(iN).b)), ...
        'VariableNames',{'Caso','MaxDiff_alpha_deg','MaxDiff_beta_deg','RMSDiff_alpha_deg','RMSDiff_beta_deg'});
    Tn = [Tn; row]; 
end
if ~isempty(Tn)
    disp('=== SCARTO LINEARE vs NON LINEARE ===');
    disp(Tn);
end

%% ========================================================================
%  ANALISI MATEMATICA (modello lineare nominale)
% ========================================================================
P  = ss(P_ext.A, P_ext.B, P_ext.C, P_ext.D);
Vh = chol(0.5*(sensor.R + sensor.R.'), 'lower');

H2L = load('H2_LQG_controller.mat',  'K_H2_lqg', 'P_H2');
H2S = load('H2_soft_controller.mat', 'K_H2_soft', 'Bd');
LQ  = load('LQG_Controllers.mat',    'K_lqg_2dof', 'K_lqg_int');
PH2 = H2L.P_H2;
Bd  = H2S.Bd;                                   % ingresso delle coppie sugli stati (nx x 2)

ctrl = struct('name',{'LQG 2DOF','H2-LQG','H2-soft','LQG integratore'}, ...
              'K',   {ss(LQ.K_lqg_2dof(1:2,:)), ss(H2L.K_H2_lqg), ...
                      ss(H2S.K_H2_soft),        ss(LQ.K_lqg_int(1:2,:))});
Tm = struct2table([res{:}]);
Tm = movevars(Tm, 'Caso', 'Before', 1);
Jfin = Tm.J_H2(isfinite(Tm.J_H2));
if ~isempty(Jfin), Tm.J_su_min = Tm.J_H2 / min(Jfin); end

disp('=== ANALISI MATEMATICA (modello lineare nominale) ===');
disp(Tm);
disp(['Legenda: Ms/Mt = picchi di S e T all''ingresso dell''impianto; DiskMargin = margine a disco ', ...
      '(piu'' alto = meglio); wB_S = banda di S a -3 dB [rad/s]; wB_T = banda di T (sigma minimo) a -3 dB; ', ...
      'Mp_T = picco del trasferimento r -> angoli; KS_inf = picco rumore IMU -> comando; ', ...
      'RMSu_rumore = RMS analitico del comando per rumore bianco (solo confronto relativo, Inf = feedthrough); ', ...
      'Off_* = offset atteso a regime per i disturbi di progetto (analitico); Tau_lento = costante di tempo del polo ', ...
      'piu'' lento; J_H2 = costo LQG/H2 sul problema comune (minimo = H2/LQG).']);

%% ========================================================================
%  SALVATAGGIO E GRAFICI (cartella nuova, niente sovrascritture)
%% ========================================================================
% --- tabelle CSV
writetable(Tt, fullfile(cfg.outDir, 'metriche_tempo.csv'));
writetable(Tm, fullfile(cfg.outDir, 'analisi_matematica.csv'));
if ~isempty(Tn), writetable(Tn, fullfile(cfg.outDir, 'scarto_lin_nl.csv')); end
if ~isempty(Te), writetable(Te, fullfile(cfg.outDir, 'equivalenza.csv')); end

% --- serie temporali ricampionate 
Tser = table(tc, 'VariableNames', {'t_s'});
for k = 1:nC
    nm = matlab.lang.makeValidName(cases(k).label);
    Tser.([nm '_alpha_deg']) = rad2deg(R(k).a);
    Tser.([nm '_beta_deg'])  = rad2deg(R(k).b);
    if ~isempty(R(k).u)
        Tser.([nm '_dF1_N']) = R(k).u(:,1);
        Tser.([nm '_dF2_N']) = R(k).u(:,2);
    end
end
writetable(Tser, fullfile(cfg.outDir, 'serie_temporali.csv'));

% --- un unico file Excel con un foglio per tabella
try
    xl = fullfile(cfg.outDir, 'risultati.xlsx');
    writetable(Tt, xl, 'Sheet', 'metriche_tempo');
    writetable(Tm, xl, 'Sheet', 'analisi_matematica');
    if ~isempty(Tn), writetable(Tn, xl, 'Sheet', 'scarto_lin_nl'); end
    if ~isempty(Te), writetable(Te, xl, 'Sheet', 'equivalenza'); end
catch ME
    warning('File Excel non scritto: %s', ME.message);
end

save(fullfile(cfg.outDir, 'risultati.mat'), 'Tt', 'Tm', 'Tn', 'Te', 'Tser', 'cfg', 'cases');

haveU = any(cellfun(@(x) ~isempty(x), {R.u}));
groups = unique({cases.group}, 'stable');

% --- risposte nel tempo per gruppo (lineare / non lineare)
for g = 1:numel(groups)
    idx = find(strcmp({cases.group}, groups{g}));
    fig = figure('Name', ['Tempo - ' groups{g}], 'Position',[100 100 1100 700]);
    nr = 1 + haveU;
    subplot(nr,2,1); hold on; grid on; title('\alpha [deg]');
    for k = idx, plot(tc, rad2deg(R(k).a), 'LineWidth',1.2); end
    plot(tc, rad2deg(cfg.refAlpha)*(tc>=cfg.tRef), 'k--'); xline(cfg.tDist, ':'); xlabel('t [s]');
    subplot(nr,2,2); hold on; grid on; title('\beta [deg]');
    for k = idx, plot(tc, rad2deg(R(k).b), 'LineWidth',1.2); end
    xline(cfg.tDist, ':'); xlabel('t [s]');
    if haveU
        subplot(nr,2,3); hold on; grid on; title('\DeltaF_1 [N]');
        for k = idx, if ~isempty(R(k).u), plot(tc, R(k).u(:,1)); end, end, xlabel('t [s]');
        subplot(nr,2,4); hold on; grid on; title('\DeltaF_2 [N]');
        for k = idx, if ~isempty(R(k).u), plot(tc, R(k).u(:,2)); end, end, xlabel('t [s]');
    end
    subplot(nr,2,1); legend({cases(idx).label}, 'Location','best');
    sgtitle(groups{g});
    saveFig(fig, cfg.outDir, ['tempo_' groups{g}]);

    % zoom sul disturbo
    fig = figure('Name', ['Disturbo - ' groups{g}], 'Position',[100 100 1100 450]);
    subplot(1,2,1); hold on; grid on; title('\alpha dopo il disturbo [deg]');
    for k = idx, plot(tc, rad2deg(R(k).a), 'LineWidth',1.2); end
    xlim([cfg.tDist-1 cfg.tStop]); xlabel('t [s]');
    subplot(1,2,2); hold on; grid on; title('\beta dopo il disturbo [deg]');
    for k = idx, plot(tc, rad2deg(R(k).b), 'LineWidth',1.2); end
    xlim([cfg.tDist-1 cfg.tStop]); xlabel('t [s]');
    legend({cases(idx).label}, 'Location','best'); sgtitle(['Disturbo - ' groups{g}]);
    saveFig(fig, cfg.outDir, ['disturbo_' groups{g}]);
end

% --- lineare contro non lineare, controllore per controllore
if ~isempty(Tn)
    fig = figure('Name','Lineare vs non lineare', 'Position',[100 100 1300 600]);
    for d = 1:numel(defs)
        iL = find(strcmp({cases.label}, defs(d).label), 1);
        iN = find(strcmp({cases.label}, [defs(d).label ' (NL)']), 1);
        if isempty(iL) || isempty(iN), continue; end
        subplot(2,numel(defs),d); hold on; grid on; title(['\alpha - ' defs(d).label]);
        plot(tc, rad2deg(R(iL).a), tc, rad2deg(R(iN).a), '--'); xlabel('t [s]');
        subplot(2,numel(defs),numel(defs)+d); hold on; grid on; title(['\beta - ' defs(d).label]);
        plot(tc, rad2deg(R(iL).b), tc, rad2deg(R(iN).b), '--'); xlabel('t [s]');
        if d == 1, legend('lineare','non lineare','Location','best'); end
    end
    saveFig(fig, cfg.outDir, 'lineare_vs_nonlineare');
end

% --- barre delle metriche principali
fig = figure('Name','Metriche a barre', 'Position',[100 100 1200 700]);
labs = categorical({cases.label}, {cases.label});
bv   = {'PicoBeta_dist_deg','Offset_alpha_deg','IAE_dist_alpha_degs'};
bt   = {'Picco \beta dopo il disturbo [deg]', ...
        'Offset \alpha a regime [deg]','IAE \alpha dopo il disturbo [deg s]'};
for i = 1:3
    subplot(2,2,i); bar(labs, Tt.(bv{i})); grid on; title(bt{i});
end
saveFig(fig, cfg.outDir, 'metriche_barre');

% --- S e T all'ingresso dell'impianto
okS  = ~cellfun(@isempty, Sl);
idxS = find(okS);
idxS = idxS(:).';

% (grafici con assi normali, non con i "chart" di sigma/step, cosi' savefig funziona)
fig = figure('Name','S e T all''ingresso', 'Position',[100 100 1100 450]);
subplot(1,2,1); hold on; grid on;
for i = idxS
    sv = sigma(Sl{i}, w);                      
    plot(w, 20*log10(sv(1,:)), 'LineWidth',1.3);
end
set(gca,'XScale','log'); xlabel('\omega [rad/s]'); ylabel('\sigma_{max}(S) [dB]');
title('S all''ingresso'); legend({ctrl(idxS).name}, 'Location','best');
subplot(1,2,2); hold on; grid on;
for i = idxS
    sv = sigma(Tl{i}, w);
    plot(w, 20*log10(sv(1,:)), 'LineWidth',1.3);
end
set(gca,'XScale','log'); xlabel('\omega [rad/s]'); ylabel('\sigma_{max}(T) [dB]');
title('T all''ingresso');
saveFig(fig, cfg.outDir, 'matematica_S_T');

% --- risposta analitica ai disturbi di progetto (modello lineare nominale)
fig = figure('Name','Disturbo analitico', 'Position',[100 100 1100 450]);
tv = (0:0.05:cfg.tStop - cfg.tDist).';
for i = idxS
    y = step(Gdl{i}, tv);                       % N x 2: alpha, beta
    subplot(1,2,1); hold on; plot(tv, rad2deg(y(:,1)), 'LineWidth',1.3);
    subplot(1,2,2); hold on; plot(tv, rad2deg(y(:,2)), 'LineWidth',1.3);
end
subplot(1,2,1); grid on; title('\alpha'); xlabel('t dopo il disturbo [s]'); ylabel('[deg]');
legend({ctrl(idxS).name}, 'Location','best');
subplot(1,2,2); grid on; title('\beta'); xlabel('t dopo il disturbo [s]'); ylabel('[deg]');
sgtitle('Risposta analitica a d_\alpha e d_\beta contemporanei');
saveFig(fig, cfg.outDir, 'matematica_disturbo_analitico');

fprintf('\nTabelle e figure salvate in %s:\n', cfg.outDir);
dd = dir(cfg.outDir);
fprintf('  %s\n', dd(~[dd.isdir]).name);
diary off

%% =========================== FUNZIONI LOCALI =============================
function tf = hasVar(name)
    tf = evalin('base', sprintf('exist(''%s'',''var'')', name)) == 1;
    if ~tf
        tf = ~isempty(findInOuts(name));
    end
end

function v = findInOuts(name)
    v = [];
    w = evalin('base','whos');
    for i = 1:numel(w)
        if strcmp(w(i).class, 'Simulink.SimulationOutput')
            try
                o = evalin('base', w(i).name);
                if any(strcmp(o.getElementNames, name))
                    v = o.get(name);
                    return;
                end
            catch
            end
        end
    end
end

function v = fetch(name)
    if evalin('base', sprintf('exist(''%s'',''var'')', name)) == 1
        v = evalin('base', name);
        return;
    end
    v = findInOuts(name);
    assert(~isempty(v), 'Segnale "%s" non trovato nel workspace.', name);
end

function [t,y] = getsig(name)
    s = fetch(name);
    if isa(s,'timeseries')
        t = s.Time;  y = squeeze(s.Data);
    elseif isstruct(s) && isfield(s,'time')
        t = s.time;  y = squeeze(s.signals.values);
    else
        t = fetch('tout');  y = squeeze(s);
    end
    t = t(:);
    if isvector(y), y = y(:); end
    if size(y,1) ~= numel(t), y = y.'; end
end

function saveFig(fig, folder, name)
    exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', 200);
    try
        exportgraphics(fig, fullfile(folder, [name '.pdf']), 'ContentType', 'vector');
    catch
    end
    try
        savefig(fig, fullfile(folder, [name '.fig']));
    catch
        try
            saveas(fig, fullfile(folder, [name '.fig']));
        catch
            warning('Figura "%s": formato .fig non salvato (PNG e PDF si'').', name);
        end
    end
end

function [s, Su, Tu, Gdd] = analyzeCtrl(Kfull, P, PH2, Vh, Bd, dAmp, w)
% Kfull: 2x5, ingressi [r(2); y_imu(3)], uscite u(2); u = Kr*r + Ky*y
    Kfull = ss(Kfull.A, Kfull.B, Kfull.C, Kfull.D);          % toglie i nomi dei canali
    Kr = Kfull(:,1:2);
    Ky = Kfull(:,3:5);
    nu = 2;  ny = size(P.C,1);  nx = size(P.A,1);
    I2 = ss(eye(2));

    s = struct('Ms',NaN,'Mt',NaN,'DiskMargin',NaN,'wB_S',NaN,'wB_T',NaN,'Mp_T',NaN, ...
               'KS_inf',NaN,'RMSu_rumore',NaN,'Off_alpha_deg',NaN,'Off_beta_deg',NaN, ...
               'Tau_lento_s',NaN,'J_H2',NaN);
    Su = [];  Tu = [];  Gdd = [];

    % anello all'ingresso dell'impianto
    Lu = -Ky*P;
    S0 = feedback(I2, Lu);
    if ~isstable(S0), warning('Anello chiuso instabile per un controllore.'); return; end
    Su = S0;  Tu = feedback(Lu, I2);
    s.Ms = norm(Su, inf);
    s.Mt = norm(Tu, inf);
    try
        DM = diskmargin(Lu);                    % margine a disco multi-anello
        s.DiskMargin = DM.DiskMargin;
    catch
    end
    sv = sigma(Su, w);
    i3 = find(sv(1,:) > 1/sqrt(2), 1, 'first');
    if ~isempty(i3), s.wB_S = w(i3); end

    % rumore IMU -> comando
    Tvu = feedback(Ky, P, +1);
    s.KS_inf = norm(Tvu, inf);
    if norm(Tvu.D) > 1e-9, s.RMSu_rumore = Inf; else, s.RMSu_rumore = norm(Tvu*Vh, 2); end

    % anello chiuso con ingressi [v (somma a u); d (coppie)] e uscite [alpha beta imu]
    iA = find(any(abs(P.C) > 1e-12, 1));
    C_ang = zeros(2,nx);  C_ang(1,iA(1)) = 1;  C_ang(2,iA(2)) = 1;
    Gd  = ss(P.A, [P.B Bd], [C_ang; P.C], zeros(2+ny,nu+2));
    Sel = [zeros(ny,2), eye(ny)];
    Gcl = feedback(Gd, Ky*Sel, 1:nu, 1:(2+ny), +1);

    % riferimento -> angoli (r passa da Kr)
    Tr = Gcl(1:2,1:nu)*Kr;
    s.Mp_T = norm(Tr, inf);
    svT = sigma(Tr, w);
    iT = find(svT(end,:) < 1/sqrt(2), 1, 'first');
    if ~isempty(iT), s.wB_T = w(iT); end

    % disturbi di coppia -> angoli
    Gdd = Gcl(1:2, nu+1:nu+2) * dAmp;            % 2x1: risposta ai disturbi di progetto
    off = dcgain(Gcl(1:2, nu+1:nu+2)) * dAmp;
    s.Off_alpha_deg = rad2deg(off(1));
    s.Off_beta_deg  = rad2deg(off(2));
    poli = pole(Gcl);
    s.Tau_lento_s = -1/max(real(poli));

    % costo LQG/H2 sul problema comune
    try
        CL = lft(PH2, Ky);
        if isstable(CL) && norm(CL.D) < 1e-9, s.J_H2 = norm(CL, 2); else, s.J_H2 = Inf; end
    catch
    end
end

function m = metrics(t, a, b, u, cfg)
    ref = cfg.refAlpha;
    iTrk   = t >= cfg.tRef & t < cfg.tDist;
    iPre   = t >= cfg.tDist-1 & t < cfg.tDist;
    iPost  = t >= cfg.tDist;
    iTail  = t >= cfg.tStop - cfg.tail;
    iReg   = t >= cfg.tRef + 5 & t < cfg.tDist;     % regime prima del disturbo

    at = a(iTrk);  tt = t(iTrk) - cfg.tRef;
    m.OS_alpha_pct     = max(0, (max(at) - ref)/abs(ref)*100);
    m.Tr_alpha_s       = riseTime(tt, at, ref);
    m.PicoBeta_trk_deg = rad2deg(max(abs(b(iTrk) - cfg.refBeta)));
    m.IAE_trk_alpha_degs = rad2deg(trapz(tt, abs(at - ref)));

    aPre = mean(a(iPre));  bPre = mean(b(iPre));
    m.PicoAlpha_dist_deg = rad2deg(max(abs(a(iPost) - aPre)));
    m.PicoBeta_dist_deg  = rad2deg(max(abs(b(iPost) - bPre)));
    m.Trientro_alpha_s   = settleTime(t(iPost) - cfg.tDist, a(iPost), aPre, cfg.bandDist);
    m.Trientro_beta_s    = settleTime(t(iPost) - cfg.tDist, b(iPost), bPre, cfg.bandDist);
    m.Offset_alpha_deg   = rad2deg(mean(a(iTail)) - aPre);
    m.Offset_beta_deg    = rad2deg(mean(b(iTail)) - bPre);
    m.IAE_dist_alpha_degs = rad2deg(trapz(t(iPost), abs(a(iPost) - aPre)));
    m.IAE_dist_beta_degs  = rad2deg(trapz(t(iPost), abs(b(iPost) - bPre)));

    m.RMS_alpha_reg_deg = rad2deg(rms(a(iReg) - mean(a(iReg))));   % fluttuazione a regime
    m.RMS_beta_reg_deg  = rad2deg(rms(b(iReg) - mean(b(iReg))));

    if isempty(u)
        m.RMS_u1_N = NaN;  m.RMS_u2_N = NaN;  m.Picco_u1_N = NaN;  m.Picco_u2_N = NaN;
        m.Saturazione_pct = NaN;
    else
        m.RMS_u1_N   = rms(u(iReg,1) - mean(u(iReg,1)));
        m.RMS_u2_N   = rms(u(iReg,2) - mean(u(iReg,2)));
        m.Picco_u1_N = max(abs(u(:,1)));
        m.Picco_u2_N = max(abs(u(:,2)));
        sat1 = u(:,1) >= cfg.limits(1,2)*0.999 | u(:,1) <= cfg.limits(1,1)*0.999;
        sat2 = u(:,2) >= cfg.limits(2,2)*0.999 | u(:,2) <= cfg.limits(2,1)*0.999;
        m.Saturazione_pct = 100*mean(sat1 | sat2);
    end
end

function Ts = settleTime(t, y, yfin, band)
    out = abs(y - yfin) > band;
    if ~any(out),      Ts = 0;
    elseif out(end),   Ts = NaN;
    else,              Ts = t(find(out,1,'last') + 1);
    end
end

function Tr = riseTime(t, y, yfin)
    i10 = find(y >= 0.1*yfin, 1);
    i90 = find(y >= 0.9*yfin, 1);
    if isempty(i10) || isempty(i90), Tr = NaN; else, Tr = t(i90) - t(i10); end
end