%% ========================================================================
%  PLOT RISULTATI SIMULAZIONE H-INFINITY
% ========================================================================

close all
clc

%% ------------------------------------------------------------------------
% ESTRAZIONE DATI LINEARI
% -------------------------------------------------------------------------

alpha_int  = out.alpha_hinf_int;
beta_int   = out.beta_hinf_int;

alpha_meas = out.alpha_meas_hinf;
beta_meas  = out.beta_meas_hinf;

u_cmd      = out.u_cmd_hinf;

%% ------------------------------------------------------------------------
% CONVERSIONE DEI SEGNALI IN VETTORI
% -------------------------------------------------------------------------
% Funzione locale che gestisce:
%   - timeseries
%   - struct con Time/Signals
%   - vettori numerici
%   - SimulationData.Dataset

[alpha_int,  t_alpha] = extractSignal(alpha_int);
[beta_int,  t_beta]   = extractSignal(beta_int);

[alpha_meas, t_am]    = extractSignal(alpha_meas);
[beta_meas,  t_bm]    = extractSignal(beta_meas);

[u_cmd, t_u]          = extractSignal(u_cmd);

%% ------------------------------------------------------------------------
% CONVERSIONE IN GRADI
% -------------------------------------------------------------------------

alpha_int_deg  = rad2deg(alpha_int);
beta_int_deg   = rad2deg(beta_int);

alpha_meas_deg = rad2deg(alpha_meas);
beta_meas_deg  = rad2deg(beta_meas);

%% ------------------------------------------------------------------------
% PLOT ANGOLI
% -------------------------------------------------------------------------

figure('Name','H-infinity - Angoli','Color','w');

subplot(2,1,1)

plot(t_alpha, alpha_int_deg, 'LineWidth', 1.5)

grid on
xlabel('Tempo [s]')
ylabel('\alpha [deg]')
title('Tracking \alpha')
legend('\alpha H_\infty','Location','best')

subplot(2,1,2)

plot(t_beta, beta_int_deg, 'LineWidth', 1.5)

grid on
xlabel('Tempo [s]')
ylabel('\beta [deg]')
title('Tracking \beta')
legend('\beta H_\infty','Location','best')


%% ------------------------------------------------------------------------
% PLOT CONFRONTO REALE / MISURATO
% -------------------------------------------------------------------------

figure('Name','Confronto reale-misurato','Color','w');

subplot(2,1,1)

plot(t_alpha, alpha_int_deg, 'LineWidth', 1.5)
hold on
plot(t_am, alpha_meas_deg, '--', 'LineWidth', 1.2)

grid on
xlabel('Tempo [s]')
ylabel('\alpha [deg]')
title('Confronto \alpha')
legend('\alpha_{int}','\alpha_{meas}','Location','best')

subplot(2,1,2)

plot(t_beta, beta_int_deg, 'LineWidth', 1.5)
hold on
plot(t_bm, beta_meas_deg, '--', 'LineWidth', 1.2)

grid on
xlabel('Tempo [s]')
ylabel('\beta [deg]')
title('Confronto \beta')
legend('\beta_{int}','\beta_{meas}','Location','best')


%% ------------------------------------------------------------------------
% PLOT AZIONI DI CONTROLLO
% -------------------------------------------------------------------------

figure('Name','H-infinity - Azioni di controllo','Color','w');

plot(t_u, u_cmd, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Comando')
title('Azioni di controllo u_{cmd}')

legend('u_{cmd}','Location','best')


%% ------------------------------------------------------------------------
% ERRORI DI MISURA
% -------------------------------------------------------------------------

% Interpolazione per avere gli stessi istanti temporali
alpha_meas_interp = interp1(t_am, alpha_meas_deg, t_alpha, ...
                           'linear','extrap');

beta_meas_interp  = interp1(t_bm, beta_meas_deg, t_beta, ...
                           'linear','extrap');

e_alpha_meas = alpha_int_deg - alpha_meas_interp;
e_beta_meas  = beta_int_deg  - beta_meas_interp;

figure('Name','Errori di misura','Color','w');

subplot(2,1,1)

plot(t_alpha, e_alpha_meas, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Errore [deg]')
title('Errore di misura \alpha')

subplot(2,1,2)

plot(t_beta, e_beta_meas, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Errore [deg]')
title('Errore di misura \beta')

%% ------------------------------------------------------------------------
% ESTRAZIONE DATI NON LINEARI
% -------------------------------------------------------------------------

alpha_int  = out.alpha_hinf_int_nl;
beta_int   = out.beta_hinf_int_nl;

alpha_meas = out.alpha_meas_hinf_nl;
beta_meas  = out.beta_meas_hinf_nl;

u_cmd      = out.u_cmd_hinf_nl;

%% ------------------------------------------------------------------------
% CONVERSIONE DEI SEGNALI IN VETTORI
% -------------------------------------------------------------------------
% Funzione locale che gestisce:
%   - timeseries
%   - struct con Time/Signals
%   - vettori numerici
%   - SimulationData.Dataset

[alpha_int,  t_alpha] = extractSignal(alpha_int);
[beta_int,  t_beta]   = extractSignal(beta_int);

[alpha_meas, t_am]    = extractSignal(alpha_meas);
[beta_meas,  t_bm]    = extractSignal(beta_meas);

[u_cmd, t_u]          = extractSignal(u_cmd);

%% ------------------------------------------------------------------------
% CONVERSIONE IN GRADI
% -------------------------------------------------------------------------

alpha_int_deg_nl  = rad2deg(alpha_int);
beta_int_deg_nl   = rad2deg(beta_int);

alpha_meas_deg = rad2deg(alpha_meas);
beta_meas_deg  = rad2deg(beta_meas);

%% ------------------------------------------------------------------------
% PLOT ANGOLI
% -------------------------------------------------------------------------

figure('Name','H-infinity - Angoli','Color','w');

subplot(2,1,1)

plot(t_alpha, alpha_int_deg, 'LineWidth', 1.5)

grid on
xlabel('Tempo [s]')
ylabel('\alpha [deg]')
title('Tracking \alpha')
legend('\alpha H_\infty','Location','best')

subplot(2,1,2)

plot(t_beta, beta_int_deg, 'LineWidth', 1.5)

grid on
xlabel('Tempo [s]')
ylabel('\beta [deg]')
title('Tracking \beta')
legend('\beta H_\infty','Location','best')


%% ------------------------------------------------------------------------
% PLOT CONFRONTO REALE / MISURATO
% -------------------------------------------------------------------------

figure('Name','Confronto reale-misurato','Color','w');

subplot(2,1,1)

plot(t_alpha, alpha_int_deg, 'LineWidth', 1.5)
hold on
plot(t_am, alpha_meas_deg, '--', 'LineWidth', 1.2)

grid on
xlabel('Tempo [s]')
ylabel('\alpha [deg]')
title('Confronto \alpha')
legend('\alpha_{int}','\alpha_{meas}','Location','best')

subplot(2,1,2)

plot(t_beta, beta_int_deg, 'LineWidth', 1.5)
hold on
plot(t_bm, beta_meas_deg, '--', 'LineWidth', 1.2)

grid on
xlabel('Tempo [s]')
ylabel('\beta [deg]')
title('Confronto \beta')
legend('\beta_{int}','\beta_{meas}','Location','best')


%% ------------------------------------------------------------------------
% PLOT AZIONI DI CONTROLLO
% -------------------------------------------------------------------------

figure('Name','H-infinity - Azioni di controllo','Color','w');

plot(t_u, u_cmd, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Comando')
title('Azioni di controllo u_{cmd}')

legend('u_{cmd}','Location','best')


%% ------------------------------------------------------------------------
% ERRORI DI MISURA
% -------------------------------------------------------------------------

% Interpolazione per avere gli stessi istanti temporali
alpha_meas_interp = interp1(t_am, alpha_meas_deg, t_alpha, ...
                           'linear','extrap');

beta_meas_interp  = interp1(t_bm, beta_meas_deg, t_beta, ...
                           'linear','extrap');

e_alpha_meas = alpha_int_deg - alpha_meas_interp;
e_beta_meas  = beta_int_deg  - beta_meas_interp;

figure('Name','Errori di misura','Color','w');

subplot(2,1,1)

plot(t_alpha, e_alpha_meas, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Errore [deg]')
title('Errore di misura \alpha')

subplot(2,1,2)

plot(t_beta, e_beta_meas, 'LineWidth', 1.4)

grid on
xlabel('Tempo [s]')
ylabel('Errore [deg]')
title('Errore di misura \beta')


%% ------------------------------------------------------------------------
% FUNZIONE LOCALE
% -------------------------------------------------------------------------

function [data,t] = extractSignal(x)

    % Timeseries
    if isa(x,'timeseries')
        data = x.Data;
        t = x.Time;
        return
    end

    % Struct To Workspace
    if isstruct(x)

        if isfield(x,'Time') && isfield(x,'Data')
            t = x.Time;
            data = x.Data;
            return
        end

        if isfield(x,'time') && isfield(x,'signals')
            t = x.time;
            data = x.signals.values;
            return
        end

    end

    % Numeric array
    if isnumeric(x)
        data = x;

        % Se non abbiamo il tempo, assumiamo indice campione
        t = (0:length(x)-1)';

        return
    end

    error('Formato del segnale non riconosciuto: %s',class(x));

end