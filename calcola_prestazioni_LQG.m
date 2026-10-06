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
%
% Nel Tipo 2 le metriche vengono calcolate esclusivamente nella
% fase precedente all'applicazione del disturbo.
%
% Rise time:
%   intervallo 10%-90%.
%
% Overshoot:
%   sovraelongazione rispetto al riferimento.
%
% Settling time:
%   primo istante in cui l'errore entra nella banda del 5%
%   e la media mobile dell'errore rimane entro tale banda.
%% ==========================================

t_dist = 10;

overshoot_alpha = NaN;

rise_time_alpha = NaN;

settling_time_alpha = NaN;


if (tipo_test == 1 || tipo_test == 2) && ...
        abs(r_alpha) > 1e-6


    % ----------------------------------------------------------
    % Parte precedente al disturbo
    % ----------------------------------------------------------

    idx_track = t_alpha < t_dist;

    t_track = t_alpha(idx_track);

    alpha_track = alpha_value(idx_track);


    % ----------------------------------------------------------
    % Filtro leggero nel caso del Tipo 2
    % ----------------------------------------------------------

    if tipo_test == 2

        alpha_track_metric = ...
            movmean(alpha_track,20);

    else

        alpha_track_metric = ...
            alpha_track;

    end


    % ----------------------------------------------------------
    % Overshoot
    % ----------------------------------------------------------

    overshoot_alpha = max( ...
        0, ...
        (max(alpha_track_metric)-r_alpha) ...
        / abs(r_alpha) * 100);


    % ----------------------------------------------------------
    % Rise time
    % ----------------------------------------------------------

    info_alpha = stepinfo( ...
        alpha_track_metric, ...
        t_track, ...
        'RiseTimeLimits',[0.1 0.9]);

    rise_time_alpha = ...
        info_alpha.RiseTime;

    % ----------------------------------------------------------
    % Settling time
    %
    % Banda del 10% rispetto al riferimento.
    % Il settling time e' il primo istante a partire dal quale
    % la risposta rimane nella banda fino all'applicazione
    % del disturbo a t = t_dist.
    % ----------------------------------------------------------
    
    settling_band = 0.05 * abs(r_alpha);
    
    err_settling = abs(alpha_track - r_alpha);
    
    T_hold = 0.5;
    N_hold = max(1, round(T_hold / mean(diff(t_track))));
    
    settling_time_alpha = NaN;
    
    for k = 1:(length(t_track)-N_hold+1)
    
        finestra = err_settling(k:k+N_hold-1);
    
        percentuale_in_banda = ...
            mean(finestra <= settling_band);
    
        if percentuale_in_banda >= 0.90
    
            settling_time_alpha = t_track(k);
            break
    
        end
    
    end
end
%% ==========================================
% 7. TEMPO DI RECUPERO
%
% Il recupero viene misurato dall'istante di applicazione del
% disturbo fino al primo istante in cui, dopo essere uscita
% dalla banda di accettazione, la risposta rientra nella banda
% e vi rimane per almeno 1 s.
%% ==========================================

T_rec_alpha = NaN;
T_rec_beta  = NaN;

peak_alpha_dist = NaN;
peak_beta_dist  = NaN;


if tipo_test == 2 || tipo_test == 3

    t_dist = 10;

    % ----------------------------------------------------------
    % SOGLIE
    % ----------------------------------------------------------

    soglia_alpha = deg2rad(0.3);

    
    soglia_beta = deg2rad(0.5);


    finestra_rec = 1.0;


    %% --------------------------------------
    % ALPHA
    %% --------------------------------------

    idx_post_alpha = find(t_alpha >= t_dist);

    t_post_alpha = t_alpha(idx_post_alpha);

    errore_alpha_post = ...
        abs(alpha_value(idx_post_alpha) - r_alpha);

    peak_alpha_dist = max(errore_alpha_post);


    % Individua il primo istante in cui la risposta
    % esce dalla banda dopo il disturbo
    idx_fuori_alpha = ...
        find(errore_alpha_post > soglia_alpha, 1, 'first');


    if ~isempty(idx_fuori_alpha)

        % Da questo punto in poi cerchiamo il primo rientro
        % stabile nella banda per almeno 1 s
        for k = idx_fuori_alpha:length(t_post_alpha)

            idx_end = find( ...
                t_post_alpha <= ...
                t_post_alpha(k) + finestra_rec, ...
                1, 'last');

            if isempty(idx_end)
                continue
            end

            % La finestra deve essere completa
            if t_post_alpha(idx_end) - ...
                    t_post_alpha(k) >= finestra_rec

                finestra = ...
                    errore_alpha_post(k:idx_end);

                if all(finestra <= soglia_alpha)

                    T_rec_alpha = ...
                        t_post_alpha(k) - t_dist;

                    break

                end
            end
        end

    else

        % Il disturbo non porta mai la risposta fuori dalla banda:
        % non è necessario alcun recupero.
        T_rec_alpha = 0;

    end


    %% --------------------------------------
    % BETA
    %% --------------------------------------

    idx_post_beta = find(t_beta >= t_dist);

    t_post_beta = t_beta(idx_post_beta);

    errore_beta_post = ...
        abs(beta_value(idx_post_beta) - r_beta);

    peak_beta_dist = max(errore_beta_post);


    % Individua il primo istante in cui la risposta
    % esce dalla banda dopo il disturbo
    idx_fuori_beta = ...
        find(errore_beta_post > soglia_beta, 1, 'first');


    if ~isempty(idx_fuori_beta)

        % Da questo punto in poi cerchiamo il primo rientro
        % stabile nella banda per almeno 1 s
        for k = idx_fuori_beta:length(t_post_beta)

            idx_end = find( ...
                t_post_beta <= ...
                t_post_beta(k) + finestra_rec, ...
                1, 'last');

            if isempty(idx_end)
                continue
            end

            % La finestra deve essere completa
            if t_post_beta(idx_end) - ...
                    t_post_beta(k) >= finestra_rec

                finestra = ...
                    errore_beta_post(k:idx_end);

                if all(finestra <= soglia_beta)

                    T_rec_beta = ...
                        t_post_beta(k) - t_dist;

                    break

                end
            end
        end

    else

        % Il disturbo non porta mai la risposta fuori dalla banda:
        % non è necessario alcun recupero.
        T_rec_beta = 0;

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

    fprintf( ...
    'Tempo recupero alpha              = %.4f s\n', ...
    T_rec_alpha);


    fprintf( ...
        'Tempo recupero beta               = %.4f s\n', ...
        T_rec_beta);



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
% 12. GRAFICI DISTURBANCE RECOVERY
%    Solo per tipo_test = 2 oppure 3
%% ==========================================

if tipo_test == 2 || tipo_test == 3

    %% ------------------------------------------
    % DISTURBANCE RECOVERY ALPHA
    %% ------------------------------------------

    figure( ...
        'Name',[nome ' - Alpha disturbance recovery'])

    plot( ...
        t_alpha, ...
        rad2deg(alpha_value), ...
        'LineWidth',1.5)

    hold on

    yline( ...
        rad2deg(r_alpha + soglia_alpha), ...
        '--')

    yline( ...
        rad2deg(r_alpha - soglia_alpha), ...
        '--')

    xline( ...
        t_dist, ...
        ':', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\alpha [deg]')

    legend( ...
        '\alpha', ...
        ['r_\alpha + ' num2str(rad2deg(soglia_alpha)) '°'], ...
        ['r_\alpha - ' num2str(rad2deg(soglia_alpha)) '°'], ...
        'Disturbance', ...
        'Location','best')

    title('Disturbance rejection on \alpha')


    %% ------------------------------------------
    % DISTURBANCE RECOVERY BETA
    %% ------------------------------------------

    figure( ...
        'Name',[nome ' - Beta disturbance recovery'])

    plot( ...
        t_beta, ...
        rad2deg(beta_value), ...
        'LineWidth',1.5)

    hold on

    yline( ...
        rad2deg(soglia_beta), ...
        '--')

    yline( ...
        -rad2deg(soglia_beta), ...
        '--')

    xline( ...
        t_dist, ...
        ':', ...
        'LineWidth',1)

    grid on

    xlabel('Time [s]')
    ylabel('\beta [deg]')

    legend( ...
        '\beta', ...
        ['+' num2str(rad2deg(soglia_beta)) '°'], ...
        ['-' num2str(rad2deg(soglia_beta)) '°'], ...
        'Disturbance', ...
        'Location','best')

    title('Disturbance rejection on \beta')

end


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