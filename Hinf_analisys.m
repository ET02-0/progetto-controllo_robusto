%% ========================================================================
% ANALISI COMPLETA H-INFINITY: TEORIA + VALIDAZIONE SIMULINK
% ========================================================================
clc; 

disp('============================================================');
disp(' AVVIO ANALISI COMPLETA H-INFINITY');
disp('============================================================');

I2 = eye(2);
omegaHinf = logspace(-2, 3, 500);
tStep = 0:0.01:15; % Orizzonte temporale per step teorico

% Lista controllori da analizzare 
controllersScaled = {
    K_hinfsyn_scaled
    K_mix_scaled
    K_pidcomp_scaled
};

controllersPhysical = {
    K_hinfsyn
    K_mix
    K_pidcomp
};
controllerNames = {'Hinfsyn', 'Mixsyn', 'Hinfstruct'};
nContr = length(controllerNames);


% Figure
fig_S  = figure('Name','Funzione di Sensibilità S'); hold on; grid on;
fig_KS = figure('Name','Sforzo di Controllo KS'); hold on; grid on;
fig_T  = figure('Name','Sensibilità Complementare T'); hold on; grid on;
fig_step_teorico = figure('Name','Step Risposta Teorica (Lineare)'); 


%% ========================================================================
% 2. ANALISI TEORICA (FREQUENZA E TEMPO)
% ========================================================================
for k = 1:nContr
    Ks = controllersScaled{k};
    Kp = controllersPhysical{k};
    name = controllerNames{k};
    
    fprintf('\n---> ANALISI CONTROLLORE: %s <---\n', upper(name));
    
    % --- LOOP NORMALIZZATO ---
    Ls = G_scaled * Ks;
    Ss = minreal(feedback(I2, Ls), 1e-7);
    Ts = minreal(feedback(Ls, I2), 1e-7);
    KSs = feedback(Ks, G_scaled);  
    
    
    % Plot Frequenziali
    figure(fig_S);  sigma(Ss, omegaHinf); 
    figure(fig_KS); sigma(KSs, omegaHinf);
    figure(fig_T);  sigma(Ts, omegaHinf);
    

    % --- CALCOLO MARGINI ROBUSTI (Diagnostica Veloce) ---
    polesNom = pole(Ts);
    if all(real(polesNom) < 0), disp(' [PASS] Nominal Stability (NS)'); else, disp(' [FAIL] Nominal Stability (NS)'); end
    
    %gammaNom = hinfnorm(minreal([WS*Ss; WU*KSs; WT*Ts], 1e-7));
    gammaNom = hinfnorm([WS*Ss; WU*KSs; WT*Ts]); 
    if gammaNom < 1, fprintf(' [PASS] Nominal Performance (NP) - Gamma = %.3f\n', gammaNom); else, fprintf(' [FAIL] Nominal Performance (NP) - Gamma = %.3f\n', gammaNom); end
    
    Lunc = G_uncertain_lumped_scaled * Ks;
    Su  = feedback(I2, Lunc);                        % senza minreal
    Tu  = feedback(Lunc, I2);
    KSu = feedback(Ks, G_uncertain_lumped_scaled);

    [stabMargin, ~] = robstab(feedback(Lunc, I2), robOptions('Display','off'));
    fprintf('\nRobust Stability:\n');
    fprintf('  Lower bound = %.4f\n',stabMargin.LowerBound);
    fprintf('  Upper bound = %.4f\n',stabMargin.UpperBound);
    if stabMargin.LowerBound > 1 
        fprintf(' [PASS] Robust Stability (RS) - Margin = %.2f\n', stabMargin.LowerBound);
    elseif stabMargin.UpperBound < 1
        fprintf(' [FAIL] Robust Stability (RS) - Margin = %.2f\n', stabMargin.LowerBound); 
    else
        fprintf('  [INCONCLUSIVE] Intervallo di bound non conclusivo.\n');
    end
    
    %CL_unc = [WS*Su; WU*KSu; WT*Tu]; 
    G = G_uncertain_lumped_scaled;
    G.InputName = {'u1','u2'};  
    G.OutputName = {'y1','y2'};

    K = ss(Ks);           
    K.InputName = {'e1','e2'};  
    K.OutputName = {'u1','u2'};
    
    W_S = ss(WS); 
    W_S.InputName = {'e1','e2'}; 
    W_S.OutputName = {'zS1','zS2'};
    
    W_U = ss(WU); 
    W_U.InputName = {'u1','u2'}; 
    W_U.OutputName = {'zU1','zU2'};
    
    W_T = ss(WT); 
    W_T.InputName = {'y1','y2'}; 
    W_T.OutputName = {'zT1','zT2'};
    sum1 = sumblk('e1 = r1 - y1');
    sum2 = sumblk('e2 = r2 - y2');

    CL_unc = connect(G, K, W_S, W_U, W_T, sum1, sum2, ...
        {'r1','r2'}, {'zS1','zS2','zU1','zU2','zT1','zT2'});
    %[wcGain,wcuRP,infoWCGain] = wcgain(CL_unc,wcOptions('Display','off'));
    wgrid = logspace(-2, 3, 200);
    [wcGain, wcuRP] = wcgain(ufrd(CL_unc, wgrid), ...
        wcOptions('Display','off','MussvOptions','a'));

    fprintf('\nRobust Performance:\n');
    fprintf('  Lower bound = %.4f\n',wcGain.LowerBound);
    fprintf('  Upper bound = %.4f\n',wcGain.UpperBound);
    
    if wcGain.UpperBound < 1
        fprintf('  [PASS] Robust performance garantita.\n');
    elseif wcGain.LowerBound > 1
        fprintf('  [FAIL] Robust performance sicuramente violata.\n');
    else
        fprintf('  [INCONCLUSIVE] Bound non conclusivi.\n');
    end

    % --- ANALISI TEMPORALE TEORICA ---
    Lp = G_nominal * Kp;
    Tp = minreal(feedback(Lp, I2), 1e-7);
    
    figure(fig_step_teorico);
    subplot(2,1,1); hold on; grid on;
    [y_pitch, t_p] = step(Tp(1,1) * scale_alpha, tStep);
    plot(t_p, rad2deg(y_pitch), 'LineWidth', 1.5);
    
    subplot(2,1,2); hold on; grid on;
    [y_yaw, t_y] = step(Tp(2,2) * scale_beta, tStep);
    plot(t_y, rad2deg(y_yaw), 'LineWidth', 1.5);

    %% ================================================================
    % REIEZIONE DISTURBO AERODINAMICO
    % ================================================================

    Gdcl = Su*Gdd;

    % Risposta al disturbo beta reale
    tDist = 0:0.01:15;

    d = zeros(length(tDist),2);

    idx = tDist >= aero.alpha.time;

    d(idx,1) = aero.alpha.amplitude;
    d(idx,2) = aero.beta.amplitude;

    [yDist,tDist] = lsim(Gdcl,d,tDist);

    figure;
    plot(tDist,rad2deg(yDist(:,1)),'LineWidth',1.2);
    hold on;
    plot(tDist,rad2deg(yDist(:,2)),'LineWidth',1.2);
    grid on;

    xlabel('Tempo [s]');
    ylabel('Deviazione angolare [deg]');
    title(['Reiezione disturbi aerodinamici - ' name]);

    legend('\Delta\alpha','\Delta\beta');


    %% ================================================================
    % RISPOSTA TEORICA AL RIFERIMENTO
    % ================================================================

    tStep = 0:0.01:15;

    r = zeros(length(tStep),2);

    idx = tStep >= ref.alpha.time;

    r(idx,1) = ref.alpha.final;
    r(idx,2) = ref.beta.final;

    [y,t] = lsim(Tp,r,tStep);

    figure;
    subplot(2,1,1);

    plot(t,rad2deg(y(:,1)),'LineWidth',1.5);
    hold on;
    yline(rad2deg(ref.alpha.final),'k--','Riferimento');
    grid on;

    ylabel('\alpha [deg]');
    title(['Risposta teorica - ' name]);

    subplot(2,1,2);

    plot(t,rad2deg(y(:,2)),'LineWidth',1.5);
    hold on;
    yline(rad2deg(ref.beta.final),'k--','Riferimento');
    grid on;

    ylabel('\beta [deg]');
    xlabel('Tempo [s]');

    
end

% Aggiunta dei limiti di specifica sui grafici frequenziali
figure(fig_S); sigma(inv(WS), omegaHinf, 'k--'); title('Sensibilità S(j\omega)'); legend([controllerNames, 'W_S^{-1}']);
figure(fig_KS); sigma(inv(WU), omegaHinf, 'k--'); title('Sensibilità Controllo KS(j\omega)'); legend([controllerNames, 'W_U^{-1}']);
figure(fig_T); sigma(inv(WT), omegaHinf, 'k--'); title('Sensibilità Complementare T(j\omega)'); legend([controllerNames, 'W_T^{-1}']);

% Finiture grafici teorici
figure(fig_step_teorico);
subplot(2,1,1); yline(rad2deg(scale_alpha), 'k--', 'Riferimento'); title('Step Pitch Teorico'); ylabel('Pitch [deg]'); legend(controllerNames);
subplot(2,1,2); yline(rad2deg(scale_beta), 'k--', 'Riferimento'); title('Step Yaw Teorico'); xlabel('Tempo [s]'); ylabel('Yaw [deg]');
