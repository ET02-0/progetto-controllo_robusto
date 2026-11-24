function ris = calcola_prestazioni_LQG( ...
    alpha, beta, tipo_test, nome, t_dist, ...
    r_alpha, r_beta, t_rif)

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

if ~ismember(tipo_test, [1 2 3])
    error('tipo_test deve essere 1, 2 oppure 3.')
end

% Riferimenti a gradino applicati a t_rif
r_alpha_vec = r_alpha * double(t_alpha >= t_rif);
r_beta_vec  = r_beta  * double(t_beta  >= t_rif);

%% ==========================================
% 3. ERRORI
%% ==========================================

err_alpha = r_alpha_vec - alpha_value;
err_beta  = r_beta_vec  - beta_value;

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
% 6. METRICHE DEL TRANSITORIO
%
% Overshoot: sovraelongazione rispetto al gradino.
% Rise time: intervallo 10%-90%.
% Settling time: assestamento al 5%.
%
% Tipo 1: risposta nominale senza disturbi.
% Tipo 2: tracking valutato prima del disturbo.
%
% Se il riferimento di un canale è nullo, le relative
% metriche di transitorio percentuali sono NaN.
%% ==========================================

overshoot_alpha  = NaN;
rise_time_alpha  = NaN;
settling_time_alpha = NaN;

overshoot_beta   = NaN;
rise_time_beta   = NaN;
settling_time_beta  = NaN;

if tipo_test == 1 || tipo_test == 2

    % Termine dell'intervallo considerato
    if tipo_test == 1
        t_limite = Inf;
    else
        t_limite = t_dist;
    end

    % Metriche del pitch
    [overshoot_alpha, ...
     rise_time_alpha, ...
     settling_time_alpha] = metriche_transitorio_LQG( ...
        t_alpha, alpha_value, r_alpha, ...
        t_rif, t_limite, tipo_test);

    % Metriche dello yaw
    [overshoot_beta, ...
     rise_time_beta, ...
     settling_time_beta] = metriche_transitorio_LQG( ...
        t_beta, beta_value, r_beta, ...
        t_rif, t_limite, tipo_test);

    fprintf('\n--- METRICHE DEL TRANSITORIO ---\n');

    fprintf('Overshoot alpha       = %.4f %%\n', ...
        overshoot_alpha);

    fprintf('Rise time alpha       = %.4f s\n', ...
        rise_time_alpha);

    fprintf('Settling time alpha   = %.4f s\n', ...
        settling_time_alpha);

    fprintf('Overshoot beta        = %.4f %%\n', ...
        overshoot_beta);

    fprintf('Rise time beta        = %.4f s\n', ...
        rise_time_beta);

    fprintf('Settling time beta    = %.4f s\n', ...
        settling_time_beta);

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

% Tempo di rientro nell'equilibrio dopo la perturbazione iniziale
T_ic_alpha = NaN;
T_ic_beta  = NaN;

peak_alpha_dist = NaN;
peak_beta_dist  = NaN;


if tipo_test == 2 || tipo_test == 3

    % ----------------------------------------------------------
    % SOGLIE
    % ----------------------------------------------------------

    soglia_alpha = deg2rad(0.3);

    
    soglia_beta = deg2rad(0.7);


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

%% ==========================================================
% METRICHE DELLA RISPOSTA ALLE CONDIZIONI INIZIALI
% Intervallo precedente al disturbo: 0 <= t < t_dist
% Solo per il Tipo 3
%% ==========================================================

%% ==========================================================
% METRICHE PRE-DISTURBO
% Riferimento nullo, intervallo 0 <= t < t_dist
%% ==========================================================

RMS_alpha_pre_dist  = NaN;
RMS_beta_pre_dist   = NaN;
Peak_alpha_pre_dist = NaN;
Peak_beta_pre_dist  = NaN;

if tipo_test == 3

    idx_pre_alpha = t_alpha < t_dist;
    idx_pre_beta  = t_beta  < t_dist;

    % ----------------------------------------------------------
    % TEMPO DI RIENTRO DOPO LE CONDIZIONI INIZIALI
    % Si considera esclusivamente l'intervallo precedente
    % al disturbo, quindi 0 <= t < t_dist.
    % ----------------------------------------------------------

    T_ic_alpha = tempo_rientro_stabile( ...
        t_alpha(idx_pre_alpha), ...
        abs(alpha_value(idx_pre_alpha) - r_alpha), ...
        soglia_alpha, ...
        finestra_rec);

    T_ic_beta = tempo_rientro_stabile( ...
        t_beta(idx_pre_beta), ...
        abs(beta_value(idx_pre_beta) - r_beta), ...
        soglia_beta, ...
        finestra_rec);

    if ~any(idx_pre_alpha) || ~any(idx_pre_beta)
        warning('Nessun campione disponibile prima del disturbo.');
    else

        err_alpha_pre = r_alpha - alpha_value(idx_pre_alpha);
        err_beta_pre  = r_beta  - beta_value(idx_pre_beta);

        RMS_alpha_pre_dist = rms(err_alpha_pre);
        RMS_beta_pre_dist  = rms(err_beta_pre);

        Peak_alpha_pre_dist = max(abs(err_alpha_pre));
        Peak_beta_pre_dist  = max(abs(err_beta_pre));

        fprintf('\n');
        fprintf('==================================================\n');
        fprintf('%s\n', nome);
        fprintf('RISPOSTA ALLE CONDIZIONI INIZIALI: PRIMA DEL DISTURBO\n');
        fprintf('Intervallo: 0 <= t < %.1f s\n', t_dist);
        fprintf('==================================================\n');

        fprintf('RMS errore alpha: %.4f gradi\n', ...
            rad2deg(RMS_alpha_pre_dist));

        fprintf('RMS errore beta:  %.4f gradi\n', ...
            rad2deg(RMS_beta_pre_dist));

        fprintf('Picco errore alpha: %.4f gradi\n', ...
            rad2deg(Peak_alpha_pre_dist));

        fprintf('Picco errore beta:  %.4f gradi\n', ...
            rad2deg(Peak_beta_pre_dist));

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
% 9b. METRICHE PRE-DISTURBO
% Intervallo: 7 <= t < t_dist
%
% Tipo 2: tracking con riferimento non nullo.
% Tipo 3: regolazione con riferimento nullo.
%
% Errore medio: misura l'eventuale offset residuo.
% RMS: misura l'errore complessivo, incluse le oscillazioni.
%% ==========================================

RMS_alpha_rumore = NaN;
RMS_beta_rumore  = NaN;

mean_error_alpha_pre_dist = NaN;
mean_error_beta_pre_dist  = NaN;

abs_mean_error_alpha_pre_dist = NaN;
abs_mean_error_beta_pre_dist  = NaN;

if tipo_test == 2 || tipo_test == 3

    idx_alpha_rumore = ...
        t_alpha >= 7 & t_alpha < t_dist;

    idx_beta_rumore = ...
        t_beta >= 7 & t_beta < t_dist;

    if any(idx_alpha_rumore) && any(idx_beta_rumore)

        % --------------------------------------------------
        % ERRORI NELLA FINESTRA PRE-DISTURBO
        % errore = riferimento - uscita
        % --------------------------------------------------

        err_alpha_pre = ...
            r_alpha - alpha_value(idx_alpha_rumore);

        err_beta_pre = ...
            r_beta - beta_value(idx_beta_rumore);

        % --------------------------------------------------
        % RMS DELL'ERRORE
        % --------------------------------------------------

        RMS_alpha_rumore = rms(err_alpha_pre);
        RMS_beta_rumore  = rms(err_beta_pre);

        % --------------------------------------------------
        % ERRORE MEDIO CON SEGNO
        % --------------------------------------------------

        mean_error_alpha_pre_dist = mean(err_alpha_pre);
        mean_error_beta_pre_dist  = mean(err_beta_pre);

        % Modulo dell'errore medio
        abs_mean_error_alpha_pre_dist = ...
            abs(mean_error_alpha_pre_dist);

        abs_mean_error_beta_pre_dist = ...
            abs(mean_error_beta_pre_dist);

        % --------------------------------------------------
        % STAMPA RISULTATI
        % --------------------------------------------------

        if tipo_test == 2
            fprintf('\n--- TRACKING PRE-DISTURBO [7,10) s ---\n');
        else
            fprintf('\n--- REGOLAZIONE PRE-DISTURBO [7,10) s ---\n');
        end

        fprintf('\nErrore medio:\n');

        fprintf('Alpha = %.6f rad (%.4f deg)\n', ...
            abs_mean_error_alpha_pre_dist, ...
            rad2deg(abs_mean_error_alpha_pre_dist));

        fprintf('Beta  = %.6f rad (%.4f deg)\n', ...
            abs_mean_error_beta_pre_dist, ...
            rad2deg(abs_mean_error_beta_pre_dist));

        fprintf('\nRMS errore:\n');

        fprintf('Alpha = %.6f rad (%.4f deg)\n', ...
            RMS_alpha_rumore, ...
            rad2deg(RMS_alpha_rumore));

        fprintf('Beta  = %.6f rad (%.4f deg)\n', ...
            RMS_beta_rumore, ...
            rad2deg(RMS_beta_rumore));

    else

        warning('Campioni insufficienti nella finestra [7,t_dist).');

    end
end


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

    fprintf('\n--- PICCO ERRORE DOPO IL DISTURBO ---\n');

    fprintf('Picco alpha dopo disturbo = %.6f rad (%.3f deg)\n', ...
        peak_alpha_dist, rad2deg(peak_alpha_dist));

    fprintf('Picco beta dopo disturbo  = %.6f rad (%.3f deg)\n', ...
        peak_beta_dist, rad2deg(peak_beta_dist));


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

    fprintf('\n--- RIENTRO DOPO LE CONDIZIONI INIZIALI ---\n');

    fprintf('Tempo rientro alpha = %.4f s\n', T_ic_alpha);
    fprintf('Tempo rientro beta  = %.4f s\n', T_ic_beta);

    fprintf('\n--- RMS ERRORE - ULTIMI 5 s ---\n')

    fprintf('RMS errore alpha = %.6f rad (%.3f deg)\n', ...
        RMS_alpha_last, rad2deg(RMS_alpha_last));

    fprintf('RMS errore beta  = %.6f rad (%.3f deg)\n', ...
        RMS_beta_last, rad2deg(RMS_beta_last));

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


stairs(t_alpha, r_alpha_vec, '--', 'LineWidth', 1.2)
stairs(t_beta,  r_beta_vec,  ':',  'LineWidth', 1.2)


grid on


xlabel('Tempo [s]')

ylabel('Angolo [rad]')


legend( ...
    '\alpha', ...
    '\beta', ...
    'r_\alpha(t)', ...
    'r_\beta(t)', ...
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

    yline(rad2deg(r_beta + soglia_beta), '--')
    yline(rad2deg(r_beta - soglia_beta), '--')

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
% 13b. ERRORE BETA
%% ==========================================

figure('Name', [nome ' - Errore beta'])

plot(t_beta, err_beta, 'LineWidth', 1.5)
hold on
yline(0, '--')

grid on
xlabel('Tempo [s]')
ylabel('e_\beta [rad]')

legend( ...
    'e_\beta', ...
    'e_\beta = 0', ...
    'Location', 'best')

title('Errore di tracking \beta')


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

% Riferimenti utilizzati nell'analisi
ris.r_alpha_vec = r_alpha_vec;
ris.r_beta_vec  = r_beta_vec;

% Metriche transitorie di entrambi i canali
ris.overshoot_beta = overshoot_beta;
ris.rise_time_beta = rise_time_beta;
ris.settling_time_beta = settling_time_beta;

% Metriche pre-disturbo: valori conservati in radianti
ris.RMS_alpha_pre_dist  = RMS_alpha_pre_dist;
ris.RMS_beta_pre_dist   = RMS_beta_pre_dist;
ris.Peak_alpha_pre_dist = Peak_alpha_pre_dist;
ris.Peak_beta_pre_dist  = Peak_beta_pre_dist;

ris.RMS_alpha_rumore = RMS_alpha_rumore;
ris.RMS_beta_rumore  = RMS_beta_rumore;

% Errore medio nella finestra pre-disturbo
ris.mean_error_alpha_pre_dist = mean_error_alpha_pre_dist;
ris.mean_error_beta_pre_dist  = mean_error_beta_pre_dist;

% Modulo dell'errore medio nella finestra pre-disturbo
ris.abs_mean_error_alpha_pre_dist = ...
    abs_mean_error_alpha_pre_dist;

ris.abs_mean_error_beta_pre_dist = ...
    abs_mean_error_beta_pre_dist;

ris.T_ic_alpha = T_ic_alpha;
ris.T_ic_beta  = T_ic_beta;

end

%% ================================================================
% FUNZIONE AUSILIARIA: METRICHE DEL TRANSITORIO LQG
%% ================================================================

function [overshoot, riseTime, settlingTime] = ...
    metriche_transitorio_LQG(t, y, r, t_rif, t_limite, tipo_test)

overshoot = NaN;
riseTime = NaN;
settlingTime = NaN;

% Un riferimento nullo non consente di definire
% overshoot e rise time come percentuali del gradino.
if abs(r) < 1e-9
    return
end

% Intervallo dopo il gradino e, nel Tipo 2,
% prima dell'applicazione del disturbo.
idx = t >= t_rif & t < t_limite;

t_track = t(idx);
y_track = y(idx);

if numel(t_track) < 2
    return
end

% Filtro leggero per le metriche con rumore
if tipo_test == 2
    y_metric = movmean(y_track, 20);
else
    y_metric = y_track;
end

% Normalizzazione nella direzione del gradino.
% Funziona anche per un riferimento negativo.
y_dir = sign(r) * y_metric;
r_abs = abs(r);

% Sovraelongazione percentuale
overshoot = max(0, ...
    (max(y_dir) - r_abs) / r_abs * 100);

% Rise time: tempo tra il primo raggiungimento del 10%
% e quello del 90% del valore finale del riferimento.
idx10 = find(y_dir >= 0.1*r_abs, 1, 'first');
idx90 = find(y_dir >= 0.9*r_abs, 1, 'first');

if ~isempty(idx10) && ~isempty(idx90) && idx90 >= idx10
    riseTime = t_track(idx90) - t_track(idx10);
end

% Errore rispetto al valore finale del gradino
err = abs(y_track - r);

% Banda di assestamento del 5%
banda = 0.05 * r_abs;

if tipo_test == 1

    % Assestamento classico:
    % ultimo campione fuori banda dopo il riferimento.
    idx_out = find(err > banda, 1, 'last');

    if isempty(idx_out)

        settlingTime = 0;

    elseif idx_out == numel(t_track)

        % La risposta è ancora fuori banda alla fine.
        settlingTime = NaN;

    else

        settlingTime = ...
            t_track(idx_out + 1) - t_rif;

    end

else

    % Tipo 2: criterio operativo in presenza di rumore.
    T_hold = 1.0;
    p_min = 0.90;

    for k = 1:numel(t_track)

        j = find( ...
            t_track >= t_track(k) + T_hold, ...
            1, 'first');

        if isempty(j)
            break
        end

        percentuale_in_banda = ...
            mean(err(k:j) <= banda);

        if percentuale_in_banda >= p_min

            settlingTime = t_track(k) - t_rif;
            break

        end

    end

end

end

function T = tempo_rientro_stabile(t, errore, soglia, T_hold)

T = NaN;

t = t(:);
errore = errore(:);

if numel(t) < 2
    return
end

% Istante iniziale della simulazione
t0 = t(1);

% Cerca il primo intervallo che resta interamente
% entro la banda per almeno T_hold secondi.
for k = 1:numel(t)

    idx_end = find(t <= t(k) + T_hold, 1, 'last');

    if isempty(idx_end)
        continue
    end

    % Verifica che la finestra sia completa
    if t(idx_end) - t(k) < T_hold
        continue
    end

    % Tutti i campioni devono restare entro la banda
    if all(errore(k:idx_end) <= soglia)

        T = t(k) - t0;
        return

    end
end

end