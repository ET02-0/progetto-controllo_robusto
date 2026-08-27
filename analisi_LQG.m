close all
clc

%% ==========================================
% ANALISI PRESTAZIONI LQG - SIMULINK
%
% tipo_LQG:
% 1 = LQG 1-DOF
% 2 = LQG 2-DOF senza integratore
% 3 = LQG 2-DOF con integratore
%
% tipo_test:
% 1 = tracking nominale, senza rumore/disturbo
% 2 = tracking + disturbo + rumore: r_alpha = 0.1 rad, r_beta = 0
% 3 = disturbo + rumore
%% ==========================================

tipo_LQG  = 3;
tipo_test = 2;


%% ==========================================
% 1. ESTRAZIONE DATI DA SIMULINK
%% ==========================================

switch tipo_LQG

    case 1

        alpha = out.alpha_LQG1;
        beta  = out.beta_LQG1;

        nome = 'LQG 1-DOF senza integratore';

    case 2

        alpha = out.alpha_LQG2;
        beta  = out.beta_LQG2;

        r_alpha = 0.1;
        r_beta  = 0;

        nome = 'LQG 2-DOF senza integratore';

    case 3

        alpha = out.alpha_LQG_int;
        beta  = out.beta_LQG_int;

        r_alpha = 0.1;
        r_beta  = 0;

        nome = 'LQG 2-DOF con integratore';

    otherwise
        error('tipo_LQG deve essere 1, 2 oppure 3.')
end


t_alpha     = alpha.Time;
alpha_value = alpha.Data;

t_beta     = beta.Time;
beta_value = beta.Data;


%% ==========================================
% 2. DEFINIZIONE RIFERIMENTI ED ERRORI
%% ==========================================

switch tipo_test

    case 1
        % Tracking nominale ideale
        r_alpha = 0.1;
        r_beta  = 0;

    case 2
        % Tracking + disturbo + rumore
        r_alpha = 0.1;
        r_beta  = 0;

    case 3
        % Disturbo + rumore attorno all'equilibrio
        r_alpha = 0;
        r_beta  = 0;

    otherwise
        error('tipo_test deve essere 1, 2 oppure 3.')
end


err_alpha = r_alpha - alpha_value;
err_beta  = r_beta  - beta_value;


%% ==========================================
% 3. METRICHE BASE
%% ==========================================

max_error_alpha = max(abs(err_alpha));
max_error_beta  = max(abs(err_beta));

RMS_alpha = rms(err_alpha);
RMS_beta  = rms(err_beta);

peak_alpha = max(abs(alpha_value));
peak_beta  = max(abs(beta_value));


%% ==========================================
% 4. ERRORE FINALE
%
% Media degli ultimi 2 s per ridurre
% l'effetto del rumore
%% ==========================================

N_alpha = length(alpha_value);
N_beta  = length(beta_value);

idx_ss_alpha = t_alpha >= (t_alpha(end) - 2);
idx_ss_beta  = t_beta  >= (t_beta(end)  - 2);

mean_alpha_final = mean(alpha_value(idx_ss_alpha));
mean_beta_final  = mean(beta_value(idx_ss_beta));

final_error_alpha = abs(mean_alpha_final - r_alpha);
final_error_beta  = abs(mean_beta_final  - r_beta);

%% ==========================================
% 5. METRICHE DI TRACKING
%
% Per il test 1:
% tracking pulito.
%
% Per il test 2:
% tracking + disturbo + rumore.
% Le metriche di tracking vengono calcolate
% prima dell'applicazione del disturbo.
%% ==========================================

t_dist = 10;

overshoot_alpha = NaN;
rise_time_alpha = NaN;
settling_time_alpha = NaN;

if (tipo_test == 1 || tipo_test == 2) && abs(r_alpha) > 1e-6

    % Parte di tracking prima del disturbo
    idx_track = t_alpha <= t_dist;

    t_track = t_alpha(idx_track);
    alpha_track = alpha_value(idx_track);

    % Se c'è rumore, filtriamo leggermente solo per
    % il calcolo delle metriche di tracking
    if tipo_test == 2
        alpha_track_metric = movmean(alpha_track,20);
    else
        alpha_track_metric = alpha_track;
    end

    % Overshoot di tracking
    overshoot_alpha = max(0,...
        (max(alpha_track_metric)-r_alpha) / abs(r_alpha) * 100);

    % Rise time e settling time
    info_alpha = stepinfo(alpha_track_metric,t_track,...
        'RiseTimeLimits',[0.1 0.9]);

    rise_time_alpha = info_alpha.RiseTime;
    settling_time_alpha = info_alpha.SettlingTime;

end

%% ==========================================
% 7. TEMPO DI RECUPERO
%
% Il disturbo avviene a t = 10 s.
% Dopo il disturbo il segnale deve rimanere
% entro la soglia per almeno 1 s.
%% ==========================================

T_rec_alpha = NaN;
T_rec_beta  = NaN;

peak_alpha_dist = NaN;
peak_beta_dist  = NaN;

if tipo_test == 2 || tipo_test == 3

    t_dist = 10;

    soglia_alpha = 0.01;
    soglia_beta  = 0.07;

    finestra_rec = 1.0;

    idx_post_alpha = find(t_alpha >= t_dist);
    idx_post_beta  = find(t_beta >= t_dist);

    % Picco dopo il disturbo
    peak_alpha_dist = max(abs(alpha_value(idx_post_alpha)));
    peak_beta_dist  = max(abs(beta_value(idx_post_beta)));

    errore_alpha_post = ...
        abs(alpha_value(idx_post_alpha) - r_alpha);

    errore_beta_post = ...
        abs(beta_value(idx_post_beta) - r_beta);

    N_finestra_alpha = max(1,...
        round(finestra_rec/mean(diff(t_alpha))));

    N_finestra_beta = max(1,...
        round(finestra_rec/mean(diff(t_beta))));


    % ---- Recupero alpha ----

    for k = 1:(length(errore_alpha_post)-N_finestra_alpha+1)

        finestra = ...
            errore_alpha_post(k:k+N_finestra_alpha-1);

        if all(finestra <= soglia_alpha)

            T_rec_alpha = ...
                t_alpha(idx_post_alpha(k)) - t_dist;

            break
        end
    end


    % ---- Recupero beta ----

    for k = 1:(length(errore_beta_post)-N_finestra_beta+1)

        finestra = ...
            errore_beta_post(k:k+N_finestra_beta-1);

        if all(finestra <= soglia_beta)

            T_rec_beta = ...
                t_beta(idx_post_beta(k)) - t_dist;

            break
        end
    end


    if isnan(T_rec_alpha)
        T_rec_alpha = NaN;
    end
    
    if isnan(T_rec_beta)
        T_rec_beta = NaN;
    end

end

%% ==========================================
% 8. STEADY STATE DISTURBO + RUMORE
%% ==========================================

mean_alpha_ss = NaN;
mean_beta_ss  = NaN;

std_alpha_ss = NaN;
std_beta_ss  = NaN;

peak_alpha_ss = NaN;
peak_beta_ss  = NaN;

if tipo_test == 2 || tipo_test == 3

    idx_alpha_ss = ...
        t_alpha >= (t_alpha(end) - 2);

    idx_beta_ss = ...
        t_beta >= (t_beta(end) - 2);

    alpha_ss = alpha_value(idx_alpha_ss);
    beta_ss  = beta_value(idx_beta_ss);

    mean_alpha_ss = mean(alpha_ss);
    mean_beta_ss  = mean(beta_ss);

    std_alpha_ss = std(alpha_ss);
    std_beta_ss  = std(beta_ss);

    peak_alpha_ss = max(abs(alpha_ss));
    peak_beta_ss  = max(abs(beta_ss));

end


%% ==========================================
% 9. RISULTATI PRINCIPALI
%% ==========================================

fprintf('\n==============================================\n')
fprintf(' %s\n', nome)
fprintf('==============================================\n')

fprintf('\n--- PRESTAZIONI GENERALI ---\n')

fprintf('Picco |alpha|          = %.6f rad\n',peak_alpha);
fprintf('Picco |beta|           = %.6f rad\n',peak_beta);

fprintf('Errore massimo alpha   = %.6f rad\n',max_error_alpha);
fprintf('Errore massimo beta    = %.6f rad\n',max_error_beta);

fprintf('RMS errore alpha       = %.6f rad\n',RMS_alpha);
fprintf('RMS errore beta        = %.6f rad\n',RMS_beta);

fprintf('Errore finale alpha    = %.6f rad\n',final_error_alpha);
fprintf('Errore finale beta     = %.6f rad\n',final_error_beta);


if tipo_test == 2

    fprintf('\n--- TRACKING + DISTURBO + RUMORE ---\n')

    fprintf('Overshoot alpha        = %.4f %%\n', ...
        overshoot_alpha);

    fprintf('Rise time alpha        = %.4f s\n', ...
        rise_time_alpha);

    fprintf('Settling time alpha    = %.4f s\n', ...
        settling_time_alpha);

    fprintf('Picco alpha dopo disturbo = %.6f rad\n', ...
        peak_alpha_dist);

    fprintf('Picco beta dopo disturbo  = %.6f rad\n', ...
        peak_beta_dist);

    if isnan(T_rec_alpha)
        fprintf('Tempo recupero alfa    = NON RAGGIUNTO\n');
    else
        fprintf('Tempo recupero alfa    = %.4f s\n', T_rec_alpha);
    end

    if isnan(T_rec_beta)
        fprintf('Tempo recupero beta    = NON RAGGIUNTO\n');
    else
        fprintf('Tempo recupero beta    = %.4f s\n', T_rec_beta);
    end

    fprintf('Media alpha SS         = %.6f rad\n', ...
        mean_alpha_ss);

    fprintf('Media beta SS          = %.6f rad\n', ...
        mean_beta_ss);

    fprintf('STD alpha SS           = %.6f rad\n', ...
        std_alpha_ss);

    fprintf('STD beta SS            = %.6f rad\n', ...
        std_beta_ss);

    fprintf('Picco alpha SS         = %.6f rad\n', ...
        peak_alpha_ss);

    fprintf('Picco beta SS          = %.6f rad\n', ...
        peak_beta_ss);

elseif tipo_test == 3

    fprintf('\n--- REIEZIONE DISTURBO + RUMORE ---\n')

    fprintf('Tempo recupero alpha   = %.4f s\n', ...
        T_rec_alpha);

    fprintf('Tempo recupero beta    = %.4f s\n', ...
        T_rec_beta);

    fprintf('Media alpha SS         = %.6f rad\n', ...
        mean_alpha_ss);

    fprintf('Media beta SS          = %.6f rad\n', ...
        mean_beta_ss);

    fprintf('STD alpha SS           = %.6f rad\n', ...
        std_alpha_ss);

    fprintf('STD beta SS            = %.6f rad\n', ...
        std_beta_ss);

    fprintf('Picco alpha SS         = %.6f rad\n', ...
        peak_alpha_ss);

    fprintf('Picco beta SS          = %.6f rad\n', ...
        peak_beta_ss);

end

%% ==========================================
% 12. GRAFICO TRACKING / RISPOSTA
%% ==========================================

figure('Name',[nome ' - Tracking'])

plot(t_alpha, alpha_value, 'LineWidth', 1.5)
hold on

plot(t_beta, beta_value, 'LineWidth', 1.5)

yline(r_alpha, '--')
yline(r_beta, ':')

grid on
xlabel('Tempo [s]')
ylabel('Angolo [rad]')

legend('\alpha', '\beta', ...
       'r_\alpha', 'r_\beta', ...
       'Location', 'best')

title(nome)


%% ==========================================
% 13. FILTRO DI KALMAN
%% ==========================================

% Stati reali
x_real_ts = out.xhat_out{1}.Values;

% Stati stimati
x_hat_ts = out.xhat_out{2}.Values;

% Dati
x_real = x_real_ts.Data;
x_hat  = x_hat_ts.Data;

% Tempi associati ai dati
t_real = x_real_ts.Time;
t_hat  = x_hat_ts.Time;

% Controllo che i tempi coincidano
if length(t_real) ~= length(t_hat)
    error('Stati reali e stimati hanno un numero diverso di campioni.')
end

if any(abs(t_real - t_hat) > 1e-10)
    warning('I tempi degli stati reali e stimati non coincidono perfettamente.')
end

% Per il confronto usiamo il tempo della stima
t_kalman = t_hat;
% Stati elicottero
% [alpha alpha_dot beta beta_dot]
x_hat_plant = x_hat(:,3:6);


alpha_real = x_real(:,1);
beta_real  = x_real(:,3);

alpha_hat = x_hat_plant(:,1);
beta_hat  = x_hat_plant(:,3);


err_K_alpha = alpha_real - alpha_hat;
err_K_beta  = beta_real  - beta_hat;


RMS_K_alpha = rms(err_K_alpha);
RMS_K_beta  = rms(err_K_beta);

MAX_K_alpha = max(abs(err_K_alpha));
MAX_K_beta  = max(abs(err_K_beta));


fprintf('\n--- FILTRO DI KALMAN ---\n')

fprintf('RMS errore stima alpha = %.6f rad\n', ...
    RMS_K_alpha);

fprintf('RMS errore stima beta  = %.6f rad\n', ...
    RMS_K_beta);

fprintf('MAX errore stima alpha = %.6f rad\n', ...
    MAX_K_alpha);

fprintf('MAX errore stima beta  = %.6f rad\n', ...
    MAX_K_beta);


%% ==========================================
% 14. GRAFICO MISURA vs STIMA KALMAN
%% ==========================================

alpha_meas = out.alpha_meas_LQG;
beta_meas  = out.beta_meas_LQG;

t_meas_alpha = alpha_meas.Time;
t_meas_beta  = beta_meas.Time;

alpha_meas_value = alpha_meas.Data;
beta_meas_value  = beta_meas.Data;


figure('Name','Filtro Kalman')

subplot(2,1,1)

plot(t_meas_alpha, alpha_meas_value)
hold on
plot(t_kalman, alpha_hat, 'LineWidth', 1.5)

grid on

xlabel('Tempo [s]')
ylabel('\alpha [rad]')

legend('Misura rumorosa', ...
       'Stima Kalman', ...
       'Location','best')

title('Filtro Kalman - \alpha')


subplot(2,1,2)

plot(t_meas_beta, beta_meas_value)
hold on
plot(t_kalman, beta_hat, 'LineWidth', 1.5)

grid on

xlabel('Tempo [s]')
ylabel('\beta [rad]')

legend('Misura rumorosa', ...
       'Stima Kalman', ...
       'Location','best')

title('Filtro Kalman - \beta')


%% ==========================================
% 15. GRAFICO ERRORE DI STIMA
%% ==========================================

figure('Name','Errore filtro Kalman')

subplot(2,1,1)

plot(t_kalman, err_K_alpha, ...
     'LineWidth', 1.5)

grid on

xlabel('Tempo [s]')
ylabel('e_\alpha [rad]')

title('Errore di stima Kalman - \alpha')


subplot(2,1,2)

plot(t_kalman, err_K_beta, ...
     'LineWidth', 1.5)

grid on

xlabel('Tempo [s]')
ylabel('e_\beta [rad]')

title('Errore di stima Kalman - \beta')