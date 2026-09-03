function ris = calcola_prestazioni_LQG( ...
    alpha, ...
    beta, ...
    tipo_test, ...
    nome)

%% ==========================================
% FUNZIONE ANALISI PRESTAZIONI LQG
%% ==========================================


%% ==========================================
% 1. ESTRAZIONE DATI
%% ==========================================

t_alpha     = alpha.Time(:);
alpha_value = alpha.Data(:);

t_beta     = beta.Time(:);
beta_value = beta.Data(:);


%% ==========================================
% 2. RIFERIMENTI
%% ==========================================

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


%% ==========================================
% 3. ERRORI
%% ==========================================

err_alpha = r_alpha - alpha_value;

err_beta = r_beta - beta_value;


%% ==========================================
% 4. METRICHE BASE
%% ==========================================

max_error_alpha = max(abs(err_alpha));

max_error_beta = max(abs(err_beta));


RMS_alpha = rms(err_alpha);

RMS_beta = rms(err_beta);


peak_alpha = max(abs(alpha_value));

peak_beta = max(abs(beta_value));


%% ==========================================
% 5. ERRORE FINALE
%
% Media ultimi 2 secondi
%% ==========================================

idx_ss_alpha = ...
    t_alpha >= (t_alpha(end) - 2);

idx_ss_beta = ...
    t_beta >= (t_beta(end) - 2);


mean_alpha_final = ...
    mean(alpha_value(idx_ss_alpha));

mean_beta_final = ...
    mean(beta_value(idx_ss_beta));


final_error_alpha = ...
    abs(mean_alpha_final - r_alpha);

final_error_beta = ...
    abs(mean_beta_final - r_beta);


%% ==========================================
% 6. METRICHE DI TRACKING
%% ==========================================

t_dist = 10;


overshoot_alpha = NaN;

rise_time_alpha = NaN;

settling_time_alpha = NaN;


if (tipo_test == 1 || tipo_test == 2) && ...
        abs(r_alpha) > 1e-6


    % Parte precedente al disturbo

    idx_track = t_alpha <= t_dist;


    t_track = t_alpha(idx_track);

    alpha_track = alpha_value(idx_track);


    % Filtro leggero se presente rumore

    if tipo_test == 2

        alpha_track_metric = ...
            movmean(alpha_track,20);

    else

        alpha_track_metric = ...
            alpha_track;

    end


    % Overshoot

    overshoot_alpha = max( ...
        0, ...
        (max(alpha_track_metric)-r_alpha) ...
        / abs(r_alpha) * 100);


    % Rise / settling

    info_alpha = stepinfo( ...
        alpha_track_metric, ...
        t_track, ...
        'RiseTimeLimits',[0.1 0.9]);


    rise_time_alpha = ...
        info_alpha.RiseTime;

    settling_time_alpha = ...
        info_alpha.SettlingTime;

end


%% ==========================================
% 7. TEMPO DI RECUPERO
%% ==========================================

T_rec_alpha = NaN;

T_rec_beta = NaN;


peak_alpha_dist = NaN;

peak_beta_dist = NaN;


if tipo_test == 2 || tipo_test == 3


    t_dist = 10;


    % Soglia alpha

    soglia_alpha = deg2rad(0.3);


    % Soglia beta

    if tipo_test == 2

        soglia_beta = deg2rad(1.8);

    else

        soglia_beta = deg2rad(0.5);

    end


    finestra_rec = 1.0;


    %% --------------------------------------
    % ALPHA
    %% --------------------------------------

    idx_post_alpha = ...
        find(t_alpha >= t_dist);


    errore_alpha_post = ...
        abs(alpha_value(idx_post_alpha) - r_alpha);


    peak_alpha_dist = ...
        max(errore_alpha_post);


    [peak_alpha_error, iPeakAlphaRel] = ...
        max(errore_alpha_post);


    if peak_alpha_error <= soglia_alpha

        T_rec_alpha = 0;

    else

        error_alpha_from_peak = ...
            errore_alpha_post(iPeakAlphaRel:end);


        t_alpha_from_peak = ...
            t_alpha( ...
            idx_post_alpha(iPeakAlphaRel:end));


        Nwin_alpha = max( ...
            1, ...
            round( ...
            finestra_rec / mean(diff(t_alpha))));


        for k = 1:( ...
                length(error_alpha_from_peak) ...
                - Nwin_alpha + 1)


            finestra = ...
                error_alpha_from_peak( ...
                k:k+Nwin_alpha-1);


            if all(finestra <= soglia_alpha)

                T_rec_alpha = ...
                    t_alpha_from_peak(k) - t_dist;

                break

            end

        end

    end


    %% --------------------------------------
    % BETA
    %% --------------------------------------

    idx_post_beta = ...
        find(t_beta >= t_dist);


    errore_beta_post = ...
        abs(beta_value(idx_post_beta) - r_beta);


    peak_beta_dist = ...
        max(errore_beta_post);


    [peak_beta_error, iPeakBetaRel] = ...
        max(errore_beta_post);


    if peak_beta_error <= soglia_beta

        T_rec_beta = 0;

    else

        error_beta_from_peak = ...
            errore_beta_post(iPeakBetaRel:end);


        t_beta_from_peak = ...
            t_beta( ...
            idx_post_beta(iPeakBetaRel:end));


        Nwin_beta = max( ...
            1, ...
            round( ...
            finestra_rec / mean(diff(t_beta))));


        for k = 1:( ...
                length(error_beta_from_peak) ...
                - Nwin_beta + 1)


            finestra = ...
                error_beta_from_peak( ...
                k:k+Nwin_beta-1);


            if all(finestra <= soglia_beta)

                T_rec_beta = ...
                    t_beta_from_peak(k) - t_dist;

                break

            end

        end

    end

end


%% ==========================================
% 8. STEADY STATE
%% ==========================================

mean_alpha_ss = NaN;

mean_beta_ss = NaN;


std_alpha_ss = NaN;

std_beta_ss = NaN;


peak_alpha_ss = NaN;

peak_beta_ss = NaN;


if tipo_test == 2 || tipo_test == 3


    idx_alpha_ss = ...
        t_alpha >= (t_alpha(end)-2);


    idx_beta_ss = ...
        t_beta >= (t_beta(end)-2);


    alpha_ss = ...
        alpha_value(idx_alpha_ss);


    beta_ss = ...
        beta_value(idx_beta_ss);


    mean_alpha_ss = ...
        mean(alpha_ss);


    mean_beta_ss = ...
        mean(beta_ss);


    std_alpha_ss = ...
        std(alpha_ss);


    std_beta_ss = ...
        std(beta_ss);


    peak_alpha_ss = ...
        max(abs(alpha_ss));


    peak_beta_ss = ...
        max(abs(beta_ss));

end


%% ==========================================
% 9. PRESTAZIONI ULTIMI 5 s
%% ==========================================

T_ss = 5;


idx_last_alpha = ...
    t_alpha >= (t_alpha(end)-T_ss);


idx_last_beta = ...
    t_beta >= (t_beta(end)-T_ss);


alpha_last = ...
    alpha_value(idx_last_alpha);


beta_last = ...
    beta_value(idx_last_beta);


err_alpha_last = ...
    alpha_last - r_alpha;


err_beta_last = ...
    beta_last - r_beta;


RMS_alpha_last = ...
    rms(err_alpha_last);


PeakErr_alpha_last = ...
    max(abs(err_alpha_last));


RMS_beta_last = ...
    rms(err_beta_last);


Peak_beta_last = ...
    max(abs(err_beta_last));


%% ==========================================
% 10. RISULTATI A VIDEO
%% ==========================================

fprintf('\n')
fprintf('==============================================\n')
fprintf(' %s\n',nome)
fprintf('==============================================\n')


fprintf('\n--- PRESTAZIONI GENERALI ---\n')


fprintf('Picco |alpha|          = %.6f rad\n', ...
    peak_alpha);


fprintf('Picco |beta|           = %.6f rad\n', ...
    peak_beta);


fprintf('Errore massimo alpha   = %.6f rad\n', ...
    max_error_alpha);


fprintf('Errore massimo beta    = %.6f rad\n', ...
    max_error_beta);


fprintf('RMS errore alpha       = %.6f rad\n', ...
    RMS_alpha);


fprintf('RMS errore beta        = %.6f rad\n', ...
    RMS_beta);


fprintf('Errore finale alpha    = %.6f rad\n', ...
    final_error_alpha);


fprintf('Errore finale beta     = %.6f rad\n', ...
    final_error_beta);


if tipo_test == 2


    fprintf('\n--- TRACKING + DISTURBO + RUMORE ---\n')


    fprintf('Overshoot alpha        = %.4f %%\n', ...
        overshoot_alpha);


    fprintf('Rise time alpha        = %.4f s\n', ...
        rise_time_alpha);


    fprintf('Settling time alpha    = %.4f s\n', ...
        settling_time_alpha);


    fprintf( ...
        'Picco errore alpha dopo disturbo = %.6f rad (%.3f deg)\n', ...
        peak_alpha_dist, ...
        rad2deg(peak_alpha_dist));


    fprintf( ...
        'Picco errore beta dopo disturbo  = %.6f rad (%.3f deg)\n', ...
        peak_beta_dist, ...
        rad2deg(peak_beta_dist));


    fprintf('\n--- PRESTAZIONI A REGIME - ULTIMI 5 s ---\n')


    fprintf( ...
        'RMS errore alpha = %.6f rad (%.3f deg)\n', ...
        RMS_alpha_last, ...
        rad2deg(RMS_alpha_last));


    fprintf( ...
        'Peak errore alpha = %.6f rad (%.3f deg)\n', ...
        PeakErr_alpha_last, ...
        rad2deg(PeakErr_alpha_last));


    fprintf( ...
        'RMS errore beta = %.6f rad (%.3f deg)\n', ...
        RMS_beta_last, ...
        rad2deg(RMS_beta_last));


    fprintf( ...
        'Peak errore beta = %.6f rad (%.3f deg)\n', ...
        Peak_beta_last, ...
        rad2deg(Peak_beta_last));


    fprintf( ...
        'Media alpha SS = %.6f rad\n', ...
        mean_alpha_ss);


    fprintf( ...
        'Media beta SS = %.6f rad\n', ...
        mean_beta_ss);


    fprintf( ...
        'STD alpha SS = %.6f rad\n', ...
        std_alpha_ss);


    fprintf( ...
        'STD beta SS = %.6f rad\n', ...
        std_beta_ss);


    fprintf( ...
        'Picco alpha SS = %.6f rad\n', ...
        peak_alpha_ss);


    fprintf( ...
        'Picco beta SS = %.6f rad\n', ...
        peak_beta_ss);


elseif tipo_test == 3


    fprintf('\n--- REIEZIONE DISTURBO + RUMORE ---\n')


    fprintf( ...
        'Tempo recupero alpha = %.4f s\n', ...
        T_rec_alpha);


    fprintf( ...
        'Tempo recupero beta = %.4f s\n', ...
        T_rec_beta);


    fprintf( ...
        'Media alpha SS = %.6f rad\n', ...
        mean_alpha_ss);


    fprintf( ...
        'Media beta SS = %.6f rad\n', ...
        mean_beta_ss);


    fprintf( ...
        'STD alpha SS = %.6f rad\n', ...
        std_alpha_ss);


    fprintf( ...
        'STD beta SS = %.6f rad\n', ...
        std_beta_ss);


    fprintf( ...
        'Picco alpha SS = %.6f rad\n', ...
        peak_alpha_ss);


    fprintf( ...
        'Picco beta SS = %.6f rad\n', ...
        peak_beta_ss);

end


%% ==========================================
% 11. GRAFICO TRACKING
%% ==========================================

figure( ...
    'Name',[nome ' - Tracking'])


plot( ...
    t_alpha, ...
    alpha_value, ...
    'LineWidth',1.5)

hold on


plot( ...
    t_beta, ...
    beta_value, ...
    'LineWidth',1.5)


yline(r_alpha,'--')

yline(r_beta,':')


grid on


xlabel('Tempo [s]')

ylabel('Angolo [rad]')


legend( ...
    '\alpha', ...
    '\beta', ...
    'r_\alpha', ...
    'r_\beta', ...
    'Location','best')


title(nome)

%% ==========================================
% 12b-alpha. REIEZIONE DISTURBO ALPHA
%% ==========================================

figure('Name',[nome ' - Alpha disturbance recovery'])

plot(t_alpha, alpha_value, 'LineWidth', 1.5)
hold on

yline(r_alpha + soglia_alpha, '--')
yline(r_alpha - soglia_alpha, '--')
xline(t_dist, ':')

grid on
xlabel('Tempo [s]')
ylabel('\alpha [rad]')

legend('\alpha', ...
       ['r_\alpha + ' num2str(rad2deg(soglia_alpha)) '°'], ...
       ['r_\alpha - ' num2str(rad2deg(soglia_alpha)) '°'], ...
       'Disturbo', ...
       'Location','best')

title('Reiezione del disturbo su \alpha')


%% ==========================================
% 12. RECUPERO BETA
%% ==========================================

figure( ...
    'Name',[nome ' - Beta disturbance recovery'])


plot( ...
    t_beta, ...
    beta_value, ...
    'LineWidth',1.5)

hold on


if tipo_test == 2 || tipo_test == 3

    yline( ...
        soglia_beta, ...
        '--')

    yline( ...
        -soglia_beta, ...
        '--')

    xline( ...
        t_dist, ...
        ':')

end


grid on


xlabel('Tempo [s]')

ylabel('\beta [rad]')


if tipo_test == 2 || tipo_test == 3

    legend( ...
        '\beta', ...
        ['+' num2str(rad2deg(soglia_beta)) '°'], ...
        ['-' num2str(rad2deg(soglia_beta)) '°'], ...
        'Disturbo', ...
        'Location','best')

else

    legend('\beta','Location','best')

end


title('Reiezione del disturbo su \beta')


%% ==========================================
% 13. ERRORE ALPHA
%% ==========================================

figure( ...
    'Name',[nome ' - Errore alpha'])


plot( ...
    t_alpha, ...
    err_alpha, ...
    'LineWidth',1.5)

hold on


yline(0,'--')


grid on


xlabel('Tempo [s]')

ylabel('e_\alpha [rad]')


legend( ...
    'e_\alpha', ...
    'e_\alpha = 0', ...
    'Location','best')


title('Errore di tracking \alpha')


%% ==========================================
% 14. SALVATAGGIO RISULTATI
%% ==========================================

ris.nome = nome;


ris.t_alpha = t_alpha;

ris.alpha_value = alpha_value;


ris.t_beta = t_beta;

ris.beta_value = beta_value;


ris.r_alpha = r_alpha;

ris.r_beta = r_beta;


ris.err_alpha = err_alpha;

ris.err_beta = err_beta;


ris.max_error_alpha = max_error_alpha;

ris.max_error_beta = max_error_beta;


ris.RMS_alpha = RMS_alpha;

ris.RMS_beta = RMS_beta;


ris.peak_alpha = peak_alpha;

ris.peak_beta = peak_beta;


ris.final_error_alpha = final_error_alpha;

ris.final_error_beta = final_error_beta;


ris.overshoot_alpha = overshoot_alpha;

ris.rise_time_alpha = rise_time_alpha;

ris.settling_time_alpha = settling_time_alpha;


ris.T_rec_alpha = T_rec_alpha;

ris.T_rec_beta = T_rec_beta;


ris.peak_alpha_dist = peak_alpha_dist;

ris.peak_beta_dist = peak_beta_dist;


ris.mean_alpha_ss = mean_alpha_ss;

ris.mean_beta_ss = mean_beta_ss;


ris.std_alpha_ss = std_alpha_ss;

ris.std_beta_ss = std_beta_ss;


ris.peak_alpha_ss = peak_alpha_ss;

ris.peak_beta_ss = peak_beta_ss;


ris.RMS_alpha_last = RMS_alpha_last;

ris.PeakErr_alpha_last = PeakErr_alpha_last;

ris.RMS_beta_last = RMS_beta_last;

ris.Peak_beta_last = Peak_beta_last;


end