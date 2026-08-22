%% ==========================================
% SIMULAZIONE MONTE CARLO (LQG)
% ==========================================


disp('--- Avvio Simulazione Monte Carlo per LQG ---')

N_campioni = 50;
rng('default');


max_overshoot_alpha = zeros(N_campioni,1);
ts_alpha = zeros(N_campioni,1);

max_beta = zeros(N_campioni,1);
rms_beta = zeros(N_campioni,1);
final_beta = zeros(N_campioni,1);

figure('Name','Monte Carlo LQG - Alpha e Beta');

subplot(2,1,1)
hold on
grid on

subplot(2,1,2)
hold on
grid on


for i = 1:N_campioni

    %% ==========================================
    % Estrazione campione incerto
    %% ==========================================

    Pk = usample(P_full_unc);

    [Ak,Bk,Ck,Dk] = ssdata(Pk);

    %% ==========================================
    % Aggiunta delle uscite fisiche alpha e beta
    %% ==========================================

    % Pk ha 10 stati.
    %
    % Dalla matrice C del modello abbiamo verificato:
    %
    % y_acc = -9.81 * x1
    % my    = -x3
    %
    % quindi:
    % alpha = x1
    % beta  = x3

    Ck_monitor = [
        Ck;
        1 0 0 0 0 0 0 0 0 0;   % alpha = x1
        0 0 1 0 0 0 0 0 0 0    % beta  = x3
    ];

    Dk_monitor = zeros(5,2);

    Pk_monitor = ss( ...
        Ak, ...
        Bk, ...
        Ck_monitor, ...
        Dk_monitor);

    Pk_monitor.InputName = {
        'u1'
        'u2'
    };

    Pk_monitor.OutputName = {
        'y_acc'
        'mx'
        'my'
        'alpha'
        'beta'
    };

    %% ==========================================
    % Chiusura dell'anello LQG
    %% ==========================================

    CLk_monitor = connect( ...
        Pk_monitor, ...
        K_lqg_int, ...
        {'r_alpha','r_beta'}, ...
        {'alpha','beta'});

    %% ==========================================
    % Risposta al riferimento reale del progetto
    % r_alpha = 0.1 rad
    % r_beta  = 0 rad
    %% ==========================================
    
    t = linspace(0,20,2001);
    
    r = zeros(length(t),2);
    r(:,1) = 0.1;    % riferimento alpha
    r(:,2) = 0;      % riferimento beta
    
    [y,t] = lsim(CLk_monitor,r,t);
    
    alpha_resp = squeeze(y(:,1));
    beta_resp  = squeeze(y(:,2));
    fprintf('Campione %d: max alpha = %.6f, max beta = %.12f\n', ...
    i, max(abs(alpha_resp)), max(abs(beta_resp)));

    %% ==========================================
    % Plot alpha
    %% ==========================================

    subplot(2,1,1)
    plot(t,alpha_resp,'Color',[0.7 0.7 0.7]);

    %% ==========================================
    % Plot beta
    %% ==========================================

    subplot(2,1,2)
    plot(t,beta_resp,'Color',[0.7 0.7 0.7]);

    %% Metriche alpha
    
    info_alpha = stepinfo(alpha_resp,t);
    
    max_overshoot_alpha(i) = info_alpha.Overshoot;
    ts_alpha(i) = info_alpha.SettlingTime;
    
    rms_alpha(i) = rms(alpha_resp - 0.1);
    final_alpha(i) = abs(alpha_resp(end) - 0.1);
    
    
    %% Metriche beta
    
    max_beta(i) = max(abs(beta_resp));
    rms_beta(i) = rms(beta_resp);
    final_beta(i) = abs(beta_resp(end));

end

%% ==========================================
% Risposta nominale
% ==========================================

t_nom = linspace(0,20,2001);

r_nom = zeros(length(t_nom),2);
r_nom(:,1) = 0.1;    % r_alpha
r_nom(:,2) = 0;      % r_beta

[y_nom,t_nom] = lsim(CL_monitor,r_nom,t_nom);

alpha_nom = squeeze(y_nom(:,1));
beta_nom  = squeeze(y_nom(:,2));

%% Plot nominale

subplot(2,1,1)

plot(t_nom,alpha_nom,'r','LineWidth',2)

xlabel('Tempo [s]')
ylabel('\alpha [rad]')
title('Monte Carlo - Tracking \alpha')
legend('Campioni incerti','Nominale','Location','best')


subplot(2,1,2)

plot(t_nom,beta_nom,'r','LineWidth',2)

xlabel('Tempo [s]')
ylabel('\beta [rad]')
title('Monte Carlo - Tracking \beta')
legend('Campioni incerti','Nominale','Location','best')


%% ==========================================
% Worst case
% ==========================================

% ALPHA
info_nom_alpha = stepinfo(alpha_nom,t_nom);

[os_alpha_worst,idx_os_alpha] = max(max_overshoot_alpha);
[ts_alpha_worst,idx_ts_alpha] = max(ts_alpha);


% BETA
max_beta_nom   = max(abs(beta_nom));
rms_beta_nom   = rms(beta_nom);
final_beta_nom = abs(beta_nom(end));

[rms_alpha_worst,idx_rms_alpha] = max(rms_alpha);
[final_alpha_worst,idx_final_alpha] = max(final_alpha);

[max_beta_worst,idx_max_beta] = max(max_beta);
[rms_beta_worst,idx_rms_beta] = max(rms_beta);
[final_beta_worst,idx_final_beta] = max(final_beta);


fprintf('\n============================================\n')
fprintf('        MONTE CARLO LQG - 50 CAMPIONI\n')
fprintf('============================================\n')


fprintf('\n--- ALPHA ---\n')

fprintf('Overshoot nominale : %.2f %%\n', ...
    info_nom_alpha.Overshoot);

fprintf('Overshoot worst     : %.2f %% (campione %d)\n', ...
    os_alpha_worst,idx_os_alpha);

fprintf('Ts nominale        : %.3f s\n', ...
    info_nom_alpha.SettlingTime);

fprintf('Ts worst            : %.3f s (campione %d)\n', ...
    ts_alpha_worst,idx_ts_alpha);


fprintf('\n--- BETA ---\n')

fprintf('Max |beta| nominale : %.6f rad\n', ...
    max_beta_nom);

fprintf('Max |beta| worst    : %.6f rad (campione %d)\n', ...
    max_beta_worst,idx_max_beta);

fprintf('RMS beta nominale   : %.6f rad\n', ...
    rms_beta_nom);

fprintf('RMS beta worst      : %.6f rad (campione %d)\n', ...
    rms_beta_worst,idx_rms_beta);

fprintf('|beta finale| nom.  : %.6f rad\n', ...
    final_beta_nom);

fprintf('|beta finale| worst : %.6f rad (campione %d)\n', ...
    final_beta_worst,idx_final_beta);